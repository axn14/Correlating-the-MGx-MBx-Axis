"""ML workbench: live weighted benchmarking, feature selection and algorithm refinement.

Launch (localhost only: the app serves clinical metadata, and Streamlit
otherwise listens on every network interface):
    PC:     py -3.12 -m streamlit run "D:\\KI Research\\v2\\ml_workbench\\app.py" --server.address localhost
    Server: ~/ki/.venv/bin/streamlit run ~/ki/v2/ml_workbench/app.py --server.address localhost --server.port 8501
            then open http://localhost:8501 through a tunnel: ssh -L 8501:localhost:8501 rs1000
            (VS Code Remote-SSH forwards the port automatically).
"""
import json
import os
import time

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import psutil
import streamlit as st
from joblib.externals.loky import get_reusable_executor
from scipy.stats import false_discovery_control
from sklearn.metrics import roc_curve

import engine as E

st.set_page_config(page_title='ML workbench', layout='wide')
ss = st.session_state

# ---------------------------------------------------------------------------
# Theme tokens (dataviz reference palette; first three categorical slots,
# which validate all-pairs in both modes)
# ---------------------------------------------------------------------------
PALETTE = {
    'light': dict(surface='#fcfcfb', ink='#0b0b0b', ink2='#52514e', muted='#898781',
                  grid='#e1e0d9', axis='#c3c2b7', s1='#2a78d6', s2='#eb6834', s3='#1baf7a',
                  seq=['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']),
    'dark': dict(surface='#1a1a19', ink='#ffffff', ink2='#c3c2b7', muted='#898781',
                 grid='#2c2c2a', axis='#383835', s1='#3987e5', s2='#d95926', s3='#199e70',
                 seq=['#0d366b', '#184f95', '#256abf', '#3987e5', '#6da7ec', '#9ec5f4', '#cde2fb']),
}
FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'


def pal() -> dict:
    return PALETTE['dark' if st.context.theme.type == 'dark' else 'light']


def group_color(group: str) -> str:
    P = pal()
    return {'species': P['s1'], 'polyamine_species': P['s2']}.get(group, P['s3'])


GROUP_LEGEND = {'species': 'Species', 'polyamine_species': 'Polyamine-pathway species',
                'clinical': 'Confounders'}


def styled(fig: go.Figure, title: str | None = None, height: int = 360) -> go.Figure:
    P = pal()
    fig.update_layout(
        template='none', height=height, barcornerradius=4,
        margin=dict(l=8, r=8, t=56 if title else 16, b=8),
        paper_bgcolor=P['surface'], plot_bgcolor=P['surface'],
        font=dict(family=FONT, color=P['ink2'], size=12),
        legend=dict(orientation='h', yanchor='bottom', y=1.0, x=0, font=dict(color=P['ink2']),
                    bgcolor='rgba(0,0,0,0)'),
        hoverlabel=dict(bgcolor=P['surface'], font=dict(color=P['ink'], family=FONT), bordercolor=P['axis']),
    )
    if title:
        fig.update_layout(title=dict(text=title, x=0, xanchor='left', y=0.98, yanchor='top',
                                     font=dict(color=P['ink'], size=14)))
    for ax in (fig.update_xaxes, fig.update_yaxes):
        ax(gridcolor=P['grid'], linecolor=P['axis'], zerolinecolor=P['axis'],
           tickfont=dict(color=P['muted']), title_font=dict(color=P['ink2']), automargin=True)
    return fig


def show(fig):
    st.plotly_chart(fig, width='stretch', theme=None, config={'displaylogo': False})


# ---------------------------------------------------------------------------
# Data and session state
# ---------------------------------------------------------------------------

@st.cache_resource(show_spinner='Loading v2 data and AGORA2 index...')
def data_for(holdout_frac: float, split_seed: int) -> dict:
    return E.load_workbench_data(holdout_frac, split_seed)


def init_state():
    if 'jobs' in ss:
        return
    ss.jobs, ss.runs, ss.nulls, ss.holdouts, ss.audits = [], {}, {}, [], {}
    ss.params_override, ss.focus = {}, None
    for rec in E.read_log():
        kind = rec.get('kind')
        if kind == 'run':
            ss.runs[rec['run_id']] = rec
        elif kind == 'null':
            ss.nulls.setdefault(rec['run_id'], []).extend(rec['values'])
        elif kind == 'holdout':
            ss.holdouts.append(rec)
        elif kind == 'audit':
            ss.audits[rec['key']] = rec


init_state()
N_CPU = os.cpu_count() or 2
N_PHYSICAL = psutil.cpu_count(logical=False) or N_CPU
# Cores this process may actually use (a cgroup/cpuset limit, e.g. on the
# shared server); os.cpu_count() and psutil ignore such limits.
N_ALLOWED = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else N_CPU
DEFAULT_WORKERS = max(1, min(N_PHYSICAL - 1, N_ALLOWED))


def executor(n_workers: int):
    return get_reusable_executor(max_workers=n_workers, timeout=900)


# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------

def active(kind=None, key=None):
    return [j for j in ss.jobs if not j['done'] and (kind is None or j['kind'] == kind)
            and (key is None or j['key'] == key)]


def submit(kind, key, label, futures, **extra):
    ss.jobs.append(dict(kind=kind, key=key, label=label, futures=futures, done=False,
                        harvested=set(), started=time.time(), error=None, **extra))


def submit_eval(cfg, n_workers):
    rid = E.cfg_hash(cfg)
    if rid in ss.runs:
        ss.focus = rid
        st.toast('This configuration was already evaluated; showing the cached run.')
        return
    if active('eval', rid):
        return
    ex, cj = executor(n_workers), E.canon(cfg)
    futs = [ex.submit(E.run_fold, cj, i) for i in range(E.n_folds(cfg))]
    submit('eval', rid, f"Evaluate {cfg['target'].split('_', 1)[-1]} · {E.MODEL_LABELS[cfg['model']]}",
           futs, cfg=cfg)


def submit_null(run, n_perm, n_workers, chunk=10):
    rid = run['run_id']
    if active('null', rid):
        return
    start = 10_000 + len(ss.nulls.get(rid, []))
    seeds = list(range(start, start + n_perm))
    ex, cj = executor(n_workers), E.canon(run['cfg'])
    futs = [ex.submit(E.permutation_chunk, cj, seeds[i:i + chunk]) for i in range(0, n_perm, chunk)]
    submit('null', rid, f'Permutation null ×{n_perm} · run {rid[:6]}', futs)


def submit_tuning(cfg, n_trials, n_workers):
    name = E.create_study(cfg)
    k = max(1, min(n_workers, n_trials))
    per = [n_trials // k + (1 if i < n_trials % k else 0) for i in range(k)]
    ex, cj = executor(n_workers), E.canon(cfg)
    futs = [ex.submit(E.optuna_worker, name, cj, n, i) for i, n in enumerate(per)]
    submit('tune', name, f'Tuning {n_trials} trials · {E.MODEL_LABELS[cfg["model"]]}', futs,
           seen=0, total=n_trials)
    ss.study_pick = name


def submit_audit(cfg, n_trials, n_workers):
    key = f'{E.cfg_hash(cfg)}|{n_trials}'
    if key in ss.audits or active('audit', key):
        return
    ex, cj = executor(n_workers), E.canon(cfg)
    futs = [ex.submit(E.nested_audit_fold, cj, i, n_trials) for i in range(cfg['n_splits'])]
    submit('audit', key, f'Nested-CV audit ({n_trials} trials/fold)', futs, cfg=cfg, n_trials=n_trials)


def submit_holdout(run, n_workers):
    if active('holdout', run['run_id']):
        return
    fut = executor(n_workers).submit(E.holdout_eval, E.canon(run['cfg']))
    submit('holdout', run['run_id'], f'Holdout evaluation · run {run["run_id"][:6]}', [fut], run=run)


def study_trial_count(name):
    import optuna
    try:
        st_ = optuna.load_study(study_name=name, storage=E.optuna_storage())
    except KeyError:
        return 0
    return sum(t.state.is_finished() for t in st_.trials)


def poll() -> bool:
    """Harvest finished futures. Returns True when results changed."""
    changed = False
    for j in ss.jobs:
        if j['done']:
            continue
        if j['kind'] == 'null':
            for i, f in enumerate(j['futures']):
                if i in j['harvested'] or not f.done():
                    continue
                j['harvested'].add(i)
                if f.exception() is not None:
                    j['error'] = repr(f.exception())
                    continue
                vals = f.result()
                ss.nulls.setdefault(j['key'], []).extend(vals)
                E.append_log(dict(kind='null', run_id=j['key'], values=vals, created=time.time()))
                changed = True
            j['done'] = len(j['harvested']) == len(j['futures'])
            continue
        if j['kind'] == 'tune':
            n = study_trial_count(j['key'])
            if n != j['seen']:
                j['seen'], changed = n, True
        if not all(f.done() for f in j['futures']):
            continue
        j['done'], changed = True, True
        errs = [f.exception() for f in j['futures'] if f.exception() is not None]
        if errs:
            j['error'] = repr(errs[0])
            continue
        res = [f.result() for f in j['futures']]
        if j['kind'] == 'eval':
            rec = E.assemble_run(j['cfg'], res)
            ss.runs[rec['run_id']] = rec
            ss.focus = rec['run_id']
            E.append_log(rec)
        elif j['kind'] == 'holdout':
            run = j['run']
            rec = dict(kind='holdout', run_id=run['run_id'], target=run['target'], created=time.time(),
                       holdout_frac=run['cfg']['holdout_frac'], split_seed=run['cfg']['split_seed'], **res[0])
            ss.holdouts.append(rec)
            E.append_log(rec)
        elif j['kind'] == 'audit':
            aucs = [r['auc'] for r in sorted(res, key=lambda r: r['fold'])]
            rec = dict(kind='audit', key=j['key'], run_id=E.cfg_hash(j['cfg']), n_trials=j['n_trials'],
                       fold_aucs=aucs, auc_mean=float(np.mean(aucs)), auc_sd=float(np.std(aucs)),
                       inner_best=[r['inner_best'] for r in res], created=time.time())
            ss.audits[j['key']] = rec
            E.append_log(rec)
    return changed


def job_progress(j) -> float:
    if j['kind'] == 'tune':
        return min(1.0, j['seen'] / max(1, j['total']))
    return sum(f.done() for f in j['futures']) / len(j['futures'])


def live_status():
    changed = poll()
    for j in active():
        extra = ''
        if j['kind'] == 'eval':
            done = [f.result()['auc'] for f in j['futures'] if f.done() and f.exception() is None]
            if done:
                extra = f' — running mean AUC {np.mean(done):.3f} over {len(done)} fold(s)'
        elif j['kind'] == 'null':
            extra = f" — {len(ss.nulls.get(j['key'], []))} permutations collected"
        st.progress(job_progress(j), text=f"{j['label']}{extra} · {time.time() - j['started']:.0f}s")
    for j in [j for j in ss.jobs if j['error'] and not j.get('dismissed')]:
        st.error(f"{j['label']} failed: {j['error']}")
    if changed:
        st.rerun()


# ---------------------------------------------------------------------------
# Sidebar: configuration
# ---------------------------------------------------------------------------
sb = st.sidebar
sb.title('ML workbench')
with sb.expander('Session setup', expanded=False):
    holdout_frac = st.select_slider('Locked holdout fraction', [0.15, 0.2, 0.25], value=0.2,
                                    help='Stratified by stage group, fixed across targets. '
                                         'Never used during refinement.')
    split_seed = int(st.number_input('Split seed', value=42, step=1))
    n_workers = st.slider('Worker processes', 1, N_ALLOWED, value=DEFAULT_WORKERS,
                          help='Default: physical cores − 1, capped at the cores this process may '
                               'use. Each model is single-threaded; more workers than cores slows '
                               'everything down.')

data = data_for(holdout_frac, split_seed)
targets = data['targets']
groups = data['groups']
group_n = pd.Series(groups).value_counts()


def target_label(f):
    r = targets.set_index('Feature').loc[f]
    tag = ' · polyamine' if r['is_polyamine'] else ''
    return f"{r['name']}{tag} (min q {r['min_q']:.3g})"


target = sb.selectbox('Target metabolite', targets['Feature'], format_func=target_label, key='w_target')
model = sb.selectbox('Model family', E.MODELS, format_func=E.MODEL_LABELS.get, key='w_model')

sb.subheader('Feature-group weights')
sb.caption('0 removes the group. Other values scale selection scores and act inside the model '
           '(XGBoost column sampling, LightGBM split gain, elastic-net penalty). '
           'Random forest has no native mechanism: weights act at selection only.')
gw = {}
for g, default in [('species', 1.0), ('polyamine_species', 1.0), ('clinical', 1.0)]:
    gw[g] = sb.slider(f"{E.GROUP_LABELS[g]} ({group_n.get(g, 0)})", 0.0, 3.0, default, 0.25, key=f'w_gw_{g}')
if sum(gw.values()) == 0:
    sb.error('At least one feature group needs a weight above 0.')

sb.subheader('Sample weights')
class_balance = sb.toggle('Balance classes', value=False, key='w_cb')
stage_w = {}
cols = sb.columns(2)
for i, g in enumerate(E.utils.STAGE_GROUP_ORDER):
    stage_w[g] = cols[i % 2].number_input(g.replace('_', ' '), 0.0, 5.0, 1.0, 0.25, key=f'w_sw_{g}')
if sum(stage_w.values()) == 0:
    sb.error('At least one stage group needs a weight above 0.')

sb.subheader('Feature selection (inside each training fold)')
fs_method = sb.selectbox('Method', E.FS_METHODS, format_func=E.FS_LABELS.get, key='w_fs')
fs_k, fs_params = 0, {}
if fs_method != 'none':
    fs_k = sb.slider('Keep at most k features', 5, 200, 50, 5, key='w_k')
if fs_method in ('l1', 'stability'):
    fs_params['l1_C'] = sb.select_slider('L1 strength C', [0.01, 0.03, 0.1, 0.3, 1.0], 0.1, key='w_c')
if fs_method == 'stability':
    fs_params['stability_B'] = sb.slider('Subsamples B', 10, 100, 30, 10, key='w_b')
    fs_params['stability_pi'] = sb.slider('Selection threshold π', 0.5, 0.9, 0.6, 0.05, key='w_pi')

sb.subheader('Cross-validation (development set)')
c1, c2 = sb.columns(2)
n_splits = c1.number_input('Folds', 3, 10, 5, key='w_splits')
n_repeats = c2.number_input('Repeats', 1, 5, 1, key='w_repeats')
cv_seed = int(sb.number_input('CV seed', value=42, step=1, key='w_cvseed'))

params = ss.params_override.get(model, {})
if params:
    sb.caption('Parameter overrides from tuning:')
    sb.json(params, expanded=False)
    if sb.button('Reset to defaults'):
        ss.params_override.pop(model, None)
        st.rerun()

cfg = dict(target=target, model=model, params=params, group_weights=gw, class_balance=class_balance,
           stage_weights=stage_w, fs_method=fs_method, fs_k=int(fs_k), fs_params=fs_params,
           n_splits=int(n_splits), n_repeats=int(n_repeats), cv_seed=cv_seed,
           holdout_frac=float(holdout_frac), split_seed=split_seed)
cfg_ok = sum(gw.values()) > 0 and sum(stage_w.values()) > 0
current_id = E.cfg_hash(cfg)

if sb.button('Evaluate configuration', type='primary', width='stretch', disabled=not cfg_ok):
    submit_eval(cfg, n_workers)
    st.rerun()

# ---------------------------------------------------------------------------
# Header, live status
# ---------------------------------------------------------------------------
trow = targets.set_index('Feature').loc[target]
y_all, thr = E.target_labels(data, target)
n_pos_dev = int(y_all[data['dev_pos']].sum())
st.title(trow['name'])
st.caption(
    f"Early CRC Cliff's δ {trow['early_CliffsD']:+.2f} (q {trow['early_QValue']:.3g}) · "
    f"Advanced CRC Cliff's δ {trow['advanced_CliffsD']:+.2f} (q {trow['advanced_QValue']:.3g}) · "
    f"development set {len(data['dev_pos'])} samples ({n_pos_dev} high) · locked holdout {len(data['hold_pos'])}"
    + (' · polyamine' if trow['is_polyamine'] else ''))

st.fragment(live_status, run_every=2 if active() else None)()

# Runs on the current split
runs_all = [r for r in ss.runs.values()
            if r['cfg']['holdout_frac'] == holdout_frac and r['cfg']['split_seed'] == split_seed]
if current_id in ss.runs:
    st.success(f"Current configuration evaluated: AUC {ss.runs[current_id]['auc_mean']:.3f} "
               f"(run {current_id[:6]}).")
elif active('eval', current_id):
    st.info('Current configuration is being evaluated.')
else:
    st.info('Current configuration not evaluated yet. Press **Evaluate configuration** in the sidebar.')


def run_label(r):
    return f"{E.MODEL_LABELS[r['model']]} · {r['fs_method']} · {r['run_id'][:5]}"


def perm_summary(r):
    null = ss.nulls.get(r['run_id'], [])
    if not null:
        return np.nan, 0
    return E.perm_pvalue(r['auc_mean'], null)[0], len(null)


def latest_holdout(run_id):
    hs = [h for h in ss.holdouts if h['run_id'] == run_id]
    return hs[-1] if hs else None


tabs = st.tabs(['Benchmark', 'Feature selection', 'Tuning', 'Significance', 'Holdout', 'Export'])

# ---------------------------------------------------------------------------
# Benchmark
# ---------------------------------------------------------------------------
with tabs[0]:
    st.subheader('Composite score weights')
    wc = st.columns(4)
    cw = dict(auc=wc[0].slider('AUC', 0.0, 1.0, 1.0, 0.05),
              stability=wc[1].slider('Stability (1 − fold SD)', 0.0, 1.0, 0.5, 0.05),
              parsimony=wc[2].slider('Parsimony (fewer features)', 0.0, 1.0, 0.25, 0.05),
              gain=wc[3].slider('Gain over confounder baseline', 0.0, 1.0, 1.0, 0.05))
    scope = st.radio('Leaderboard scope', ['This target', 'All targets'], horizontal=True,
                     help='AUCs from different targets are not comparable; compare within a target.')
    board_runs = [r for r in runs_all if scope == 'All targets' or r['target'] == target]

    if not board_runs:
        st.info('No runs yet for this target on this split.')
    else:
        df = pd.DataFrame(board_runs)
        scored = E.composite_scores(df, cw)
        df = df.join(scored)
        df['perm_p'], df['n_perm'] = zip(*[perm_summary(r) for r in board_runs])
        df['holdout_auc'] = [(latest_holdout(r['run_id']) or {}).get('auc', np.nan) for r in board_runs]
        df['weights'] = [
            ' '.join(f"{k[:3]}={v:g}" for k, v in r['cfg']['group_weights'].items())
            + (' · balanced' if r['cfg']['class_balance'] else '')
            + ('' if set(r['cfg']['stage_weights'].values()) == {1.0} else ' · stage-weighted')
            for r in board_runs]
        df['metabolite'] = [targets.set_index('Feature').loc[r['target'], 'name'] for r in board_runs]
        df['model_label'] = df['model'].map(E.MODEL_LABELS)
        df = df.sort_values('composite', ascending=False).reset_index(drop=True)
        view = df[['composite', 'metabolite', 'model_label', 'fs_method', 'auc_mean', 'auc_sd',
                   'baseline_auc_mean', 'gain', 'n_features_mean', 'perm_p', 'holdout_auc', 'weights', 'run_id']]
        ev = st.dataframe(
            view, hide_index=True, width='stretch', on_select='rerun', selection_mode='single-row', key='lb',
            column_config={
                'composite': st.column_config.ProgressColumn('Composite', min_value=0, max_value=1, format='%.2f'),
                'metabolite': 'Target', 'model_label': 'Model', 'fs_method': 'Selection',
                'auc_mean': st.column_config.NumberColumn('CV AUC', format='%.3f'),
                'auc_sd': st.column_config.NumberColumn('Fold SD', format='%.3f'),
                'baseline_auc_mean': st.column_config.NumberColumn('Baseline AUC', format='%.3f'),
                'gain': st.column_config.NumberColumn('Gain', format='%+.3f'),
                'n_features_mean': st.column_config.NumberColumn('Features', format='%.0f'),
                'perm_p': st.column_config.NumberColumn('Perm p', format='%.4f'),
                'holdout_auc': st.column_config.NumberColumn('Holdout AUC', format='%.3f'),
                'weights': 'Weights', 'run_id': 'Run'})
        if ev.selection.rows:
            ss.focus = df.loc[ev.selection.rows[0], 'run_id']
        st.caption('Select a row to focus it in the other tabs. Composite parts are min-max '
                   'normalised across the runs shown, so re-weighting re-ranks instantly without refitting.')

        top = df.head(12)
        P = pal()
        fig = go.Figure()
        labels = [run_label(ss.runs[rid]) for rid in top['run_id']]
        for lab, (_, r) in zip(labels, top.iterrows()):
            fig.add_trace(go.Scatter(
                x=r['fold_aucs'], y=[lab] * len(r['fold_aucs']), mode='markers', name='Model folds',
                legendgroup='m', showlegend=lab == labels[0],
                marker=dict(size=10, color=P['s1'], line=dict(width=2, color=P['surface'])),
                hovertemplate=f'{lab}<br>fold AUC %{{x:.3f}}<extra></extra>'))
            fig.add_trace(go.Scatter(
                x=r['baseline_fold_aucs'], y=[lab] * len(r['baseline_fold_aucs']), mode='markers',
                name='Confounder baseline folds', legendgroup='b', showlegend=lab == labels[0],
                marker=dict(size=9, color=P['s2'], symbol='diamond', line=dict(width=2, color=P['surface'])),
                hovertemplate=f'{lab}<br>baseline fold AUC %{{x:.3f}}<extra></extra>'))
        fig.add_vline(x=0.5, line=dict(color=P['muted'], width=1, dash='dot'),
                      annotation_text='chance', annotation_font_color=P['muted'])
        fig.update_yaxes(categoryorder='array', categoryarray=labels[::-1])
        fig.update_xaxes(title='AUC per fold')
        show(styled(fig, 'Fold AUCs, top runs by composite score', height=120 + 34 * len(labels)))

focus_run = ss.runs.get(ss.focus) if ss.focus in ss.runs else (ss.runs.get(current_id))

with tabs[0]:
    if focus_run:
        P = pal()
        y = np.asarray(focus_run['y_dev'])
        fig = go.Figure()
        for name, key, color, auc_ in [('Model', 'oof', P['s1'], focus_run['auc_mean']),
                                       ('Confounder baseline', 'baseline_oof', P['s2'],
                                        focus_run['baseline_auc_mean'])]:
            fpr, tpr, _ = roc_curve(y, focus_run[key])
            fig.add_trace(go.Scatter(x=fpr, y=tpr, mode='lines', name=f'{name} (fold-mean AUC {auc_:.3f})',
                                     line=dict(color=color, width=2),
                                     hovertemplate='FPR %{x:.2f}<br>TPR %{y:.2f}<extra>' + name + '</extra>'))
        fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode='lines', name='Chance', hoverinfo='skip',
                                 line=dict(color=P['muted'], width=1, dash='dot')))
        fig.update_xaxes(title='False-positive rate', range=[0, 1])
        fig.update_yaxes(title='True-positive rate', range=[0, 1.02], scaleanchor='x')
        show(styled(fig, f'Out-of-fold ROC · {run_label(focus_run)}', height=440))
        if focus_run['notes']:
            st.caption('Notes: ' + '; '.join(focus_run['notes']))

# ---------------------------------------------------------------------------
# Feature selection
# ---------------------------------------------------------------------------
with tabs[1]:
    if not focus_run:
        st.info('Evaluate a configuration, or select a run on the Benchmark tab.')
    elif focus_run['fs_method'] == 'none':
        st.info('The focused run uses no feature selection: every active feature is used in every fold. '
                'Pick a selection method in the sidebar and evaluate to compare.')
    else:
        P = pal()
        freq = pd.Series(focus_run['selection_freq']).sort_values(ascending=False).head(30)
        fig = go.Figure()
        order = list(freq.index)[::-1]
        for leg in ['Species', 'Polyamine-pathway species', 'Confounders']:
            feats_ = [f for f in freq.index if GROUP_LEGEND[groups[f]] == leg]
            if not feats_:
                continue
            fig.add_trace(go.Bar(
                x=freq[feats_], y=feats_, orientation='h', name=leg,
                marker=dict(color=group_color(groups[feats_[0]]), line=dict(width=0)),
                hovertemplate='%{y}<br>selected in %{x:.0%} of folds<extra>' + leg + '</extra>'))
        fig.update_yaxes(categoryorder='array', categoryarray=order, tickfont=dict(size=11))
        fig.update_xaxes(title='Share of folds that selected the feature', tickformat='.0%', range=[0, 1.05])
        fig.update_layout(bargap=0.25)
        show(styled(fig, f'Selection frequency, top 30 · {run_label(focus_run)}', height=140 + 22 * len(freq)))

    sel_runs = [r for r in runs_all if r['target'] == target and r['fs_method'] != 'none']
    if len(sel_runs) >= 1:
        P = pal()
        comp = []
        for r in sel_runs:
            by = {'Species': 0.0, 'Polyamine-pathway species': 0.0, 'Confounders': 0.0}
            for f, v in r['selection_freq'].items():
                by[GROUP_LEGEND[groups[f]]] += v
            comp.append((run_label(r), by))
        fig = go.Figure()
        for leg, key in [('Species', 'species'), ('Polyamine-pathway species', 'polyamine_species'),
                         ('Confounders', 'clinical')]:
            fig.add_trace(go.Bar(
                x=[c[1][leg] for c in comp], y=[c[0] for c in comp], orientation='h', name=leg,
                marker=dict(color=group_color(key), line=dict(width=2, color=P['surface'])),
                hovertemplate='%{y}<br>%{x:.1f} features per fold<extra>' + leg + '</extra>'))
        fig.update_layout(barmode='stack', bargap=0.3)
        fig.update_xaxes(title='Mean features selected per fold')
        show(styled(fig, 'What each selection run keeps, by feature group', height=140 + 34 * len(comp)))

    if len(sel_runs) >= 2:
        P = pal()
        sets = {run_label(r): {f for f, v in r['selection_freq'].items() if v >= 0.5} for r in sel_runs}
        labs = list(sets)
        jac = [[len(sets[a] & sets[b]) / max(1, len(sets[a] | sets[b])) for b in labs] for a in labs]
        fig = go.Figure(go.Heatmap(
            z=jac, x=labs, y=labs, zmin=0, zmax=1,
            colorscale=[[i / (len(P['seq']) - 1), c] for i, c in enumerate(P['seq'])],
            colorbar=dict(title='Jaccard', tickfont=dict(color=P['muted'])),
            hovertemplate='%{y}<br>vs %{x}<br>overlap %{z:.2f}<extra></extra>'))
        fig.update_xaxes(showgrid=False, tickangle=-30)
        fig.update_yaxes(showgrid=False, autorange='reversed')
        show(styled(fig, 'Stable-set overlap between selection runs (features kept in ≥50% of folds)',
                    height=200 + 40 * len(labs)))

# ---------------------------------------------------------------------------
# Tuning
# ---------------------------------------------------------------------------
with tabs[2]:
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    st.caption('Optuna TPE with median pruning, tuned on the development-set CV of the current '
               'configuration (target, model, feature groups, weights, selection). The best tuning score '
               'is optimistic because it picked the winner from the same folds: use the nested-CV audit '
               'or the holdout for an honest number.')
    c1, c2 = st.columns([1, 2])
    n_trials = c1.slider('Trials', 10, 200, 30, 10)
    if c1.button('Start tuning with current configuration', disabled=not cfg_ok):
        submit_tuning(cfg, n_trials, n_workers)
        st.rerun()

    storage = E.optuna_storage()
    summaries = sorted(optuna.get_all_study_summaries(storage=storage),
                       key=lambda s: s.datetime_start or pd.Timestamp(0), reverse=True)
    names = [s.study_name for s in summaries if s.study_name.startswith(target + '|')]
    if not names:
        st.info('No tuning studies for this target yet.')
    else:
        default = names.index(ss.get('study_pick')) if ss.get('study_pick') in names else 0
        pick = c2.selectbox('Study', names, index=default,
                            format_func=lambda n: f"{E.MODEL_LABELS.get(n.split('|')[1], n)} · "
                                                  f"{n.split('|')[2]} · "
                                                  f"{time.strftime('%d %b %H:%M', time.localtime(int(n.split('|')[3])))}")
        study = optuna.load_study(study_name=pick, storage=storage)
        trials = study.trials
        done = [t for t in trials if t.state == optuna.trial.TrialState.COMPLETE]
        pruned = [t for t in trials if t.state == optuna.trial.TrialState.PRUNED]
        P = pal()
        if not done:
            st.info('Waiting for the first completed trial...')
        else:
            best_so_far = np.maximum.accumulate([t.value for t in done])
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=[t.number for t in done], y=[t.value for t in done], mode='markers', name='Completed trial',
                marker=dict(size=9, color=P['s1'], line=dict(width=2, color=P['surface'])),
                hovertemplate='trial %{x}<br>CV AUC %{y:.3f}<extra></extra>'))
            if pruned:
                fig.add_trace(go.Scatter(
                    x=[t.number for t in pruned],
                    y=[max(t.intermediate_values.values()) if t.intermediate_values else None for t in pruned],
                    mode='markers', name='Pruned (last partial score)',
                    marker=dict(size=8, color='rgba(0,0,0,0)', line=dict(width=1.5, color=P['muted'])),
                    hovertemplate='trial %{x} pruned<br>partial AUC %{y:.3f}<extra></extra>'))
            fig.add_trace(go.Scatter(
                x=[t.number for t in done], y=best_so_far, mode='lines', name='Best so far',
                line=dict(color=P['s2'], width=2, shape='hv'),
                hovertemplate='best after trial %{x}: %{y:.3f}<extra></extra>'))
            fig.update_xaxes(title='Trial')
            fig.update_yaxes(title='Tuning CV AUC (optimistic)')
            show(styled(fig, f'Optimisation history · {len(done)} complete, {len(pruned)} pruned', height=380))

            m1, m2, m3 = st.columns(3)
            m1.metric('Best tuning CV AUC', f'{study.best_value:.3f}', help='Optimistic: selected on these folds.')
            m2.metric('Trials finished', f'{len(done) + len(pruned)}')
            tune_cfg = study.user_attrs.get('cfg', {})
            base_run = ss.runs.get(E.cfg_hash({**tune_cfg, 'params': {}})) if tune_cfg else None
            if base_run:
                m3.metric('Same config at defaults', f"{base_run['auc_mean']:.3f}",
                          delta=f"{study.best_value - base_run['auc_mean']:+.3f} from tuning")

            if len(done) >= 5:
                try:
                    imp = optuna.importance.get_param_importances(study)
                except (RuntimeError, ValueError):
                    imp = {}
                if imp:
                    s = pd.Series(imp).sort_values()
                    fig = go.Figure(go.Bar(x=s.values, y=s.index, orientation='h',
                                           marker=dict(color=P['s1']),
                                           hovertemplate='%{y}: %{x:.2f}<extra></extra>'))
                    fig.update_xaxes(title='Share of explained variation in CV AUC (fANOVA)')
                    show(styled(fig, 'Parameter importance', height=120 + 30 * len(s)))

                pnames = list(done[0].params)
                dims = []
                for p_ in pnames:
                    vals = [t.params[p_] for t in done]
                    log = p_ in ('learning_rate', 'reg_lambda', 'C', 'min_child_weight')
                    dims.append(dict(label=p_ + (' (log10)' if log else ''),
                                     values=np.log10(vals) if log else vals))
                dims.append(dict(label='CV AUC', values=[t.value for t in done]))
                fig = go.Figure(go.Parcoords(
                    line=dict(color=[t.value for t in done],
                              colorscale=[[i / (len(P['seq']) - 1), c] for i, c in enumerate(P['seq'])],
                              showscale=True, colorbar=dict(title='AUC', tickfont=dict(color=P['muted']))),
                    dimensions=dims, labelfont=dict(color=P['ink2']), tickfont=dict(color=P['muted']),
                    rangefont=dict(color=P['muted'])))
                show(styled(fig, 'Where the good trials sit (darker = higher AUC)', height=420))

            st.write('**Best parameters**')
            st.json(study.best_params)
            b1, b2 = st.columns(2)

            def apply_best(tcfg=tune_cfg, best=study.best_params):
                ss.params_override[tcfg['model']] = best
                ss.w_target, ss.w_model = tcfg['target'], tcfg['model']

            b1.button('Apply best parameters to the sidebar', on_click=apply_best, disabled=not tune_cfg,
                      help='Then press Evaluate to benchmark the tuned configuration.')
            audit_trials = b2.slider('Nested-CV audit: trials per outer fold', 10, 60, 20, 10)
            if b2.button('Run nested-CV audit', disabled=not tune_cfg):
                submit_audit(tune_cfg, audit_trials, n_workers)
                st.rerun()
            audits = [a for k, a in ss.audits.items() if tune_cfg and a['run_id'] == E.cfg_hash(tune_cfg)]
            for a in audits:
                st.success(f"Nested-CV audit ({a['n_trials']} trials per fold): outer AUC "
                           f"{a['auc_mean']:.3f} ± {a['auc_sd']:.3f} vs mean inner best "
                           f"{np.mean(a['inner_best']):.3f}. The outer number is the leak-free estimate "
                           "of what tuning this pipeline achieves.")

# ---------------------------------------------------------------------------
# Significance
# ---------------------------------------------------------------------------
with tabs[3]:
    st.caption('Each permutation shuffles the labels and re-runs the full pipeline (weights, in-fold '
               'selection, model at fixed parameters) on folds re-stratified for the shuffled labels. '
               'p-values are conditional on the configuration you chose and exploratory: every '
               'configuration you tried is a comparison the p-value does not account for.')
    if not focus_run:
        st.info('Evaluate a configuration, or select a run on the Benchmark tab.')
    else:
        c1, c2 = st.columns([1, 3])
        n_perm = c1.select_slider('Permutations to add', [50, 100, 200, 500, 1000], 200)
        if c1.button('Run permutations', disabled=bool(active('null', focus_run['run_id']))):
            submit_null(focus_run, n_perm, n_workers)
            st.rerun()
        null = ss.nulls.get(focus_run['run_id'], [])
        if null:
            p, (lo, hi), b, n = E.perm_pvalue(focus_run['auc_mean'], null)
            c1.metric('Permutation p', f'{p:.4f}', help=f'(b+1)/(N+1) with b={b} null AUCs ≥ observed')
            c1.caption(f'N = {n} · Monte Carlo 95% CI for the exceedance rate [{lo:.4f}, {hi:.4f}] · '
                       f'smallest possible p at this N: {1 / (n + 1):.4f}')
            P = pal()
            fig = go.Figure(go.Histogram(x=null, nbinsx=40, name='Null AUCs', marker=dict(color=P['s1']),
                                         hovertemplate='AUC %{x}<br>%{y} permutations<extra></extra>'))
            fig.add_vline(x=focus_run['auc_mean'], line=dict(color=P['s2'], width=2))
            fig.add_annotation(x=focus_run['auc_mean'], y=1, yref='paper', showarrow=False, xanchor='left',
                               text=f" observed {focus_run['auc_mean']:.3f}", font=dict(color=P['ink']))
            fig.update_xaxes(title='Mean CV AUC under shuffled labels')
            fig.update_yaxes(title='Permutations')
            fig.update_layout(bargap=0.05, showlegend=False)
            with c2:
                show(styled(fig, f'Permutation null · {run_label(focus_run)}', height=360))

    tested = [r for r in runs_all if ss.nulls.get(r['run_id'])]
    if tested:
        rows = []
        for r in tested:
            p, _, _, n = E.perm_pvalue(r['auc_mean'], ss.nulls[r['run_id']])
            rows.append(dict(Target=targets.set_index('Feature').loc[r['target'], 'name'],
                             Run=run_label(r), AUC=r['auc_mean'], Permutations=n, p=p))
        tab = pd.DataFrame(rows)
        tab['BH q'] = false_discovery_control(tab['p'].to_numpy(), method='bh')
        st.write('**All runs with a permutation null on this split** (BH across every row: each row is a test)')
        st.dataframe(tab.sort_values('p'), hide_index=True, width='stretch',
                     column_config={'AUC': st.column_config.NumberColumn(format='%.3f'),
                                    'p': st.column_config.NumberColumn(format='%.4f'),
                                    'BH q': st.column_config.NumberColumn(format='%.4f')})

# ---------------------------------------------------------------------------
# Holdout
# ---------------------------------------------------------------------------
with tabs[4]:
    peeks = [h for h in ss.holdouts if h['target'] == target
             and h['holdout_frac'] == holdout_frac and h['split_seed'] == split_seed]
    c1, c2 = st.columns([1, 2])
    c1.metric('Holdout queries for this target', len(peeks))
    if len(peeks) > 1:
        c1.warning('Each holdout result you act on makes the next one less of an unbiased estimate. '
                   'Treat anything after the first query as another validation fold.')
    if not focus_run:
        st.info('Evaluate a configuration, or select a run on the Benchmark tab.')
    else:
        c2.write(f'Focused run: **{run_label(focus_run)}** (CV AUC {focus_run["auc_mean"]:.3f})')
        ready = c2.checkbox('I have finished refining this configuration')
        if c2.button('Evaluate on the locked holdout', type='primary',
                     disabled=not ready or bool(active('holdout', focus_run['run_id']))):
            submit_holdout(focus_run, n_workers)
            st.rerun()
    if peeks:
        h = peeks[-1]
        run = ss.runs.get(h['run_id'])
        st.subheader('Latest holdout result')
        m1, m2, m3, m4 = st.columns(4)
        m1.metric('Holdout AUC', f"{h['auc']:.3f}", help='Model refit on the whole development set.')
        m2.metric('95% bootstrap CI', f"{h['ci_low']:.3f} – {h['ci_high']:.3f}")
        m3.metric('Confounder baseline', f"{h['baseline_auc']:.3f}", delta=f"{h['auc'] - h['baseline_auc']:+.3f}")
        if run:
            m4.metric('CV AUC (development)', f"{run['auc_mean']:.3f}", delta=f"{h['auc'] - run['auc_mean']:+.3f}")
        P = pal()
        fpr, tpr, _ = roc_curve(h['y'], h['proba'])
        fig = go.Figure(go.Scatter(x=fpr, y=tpr, mode='lines', name='Holdout',
                                   line=dict(color=P['s1'], width=2),
                                   hovertemplate='FPR %{x:.2f}<br>TPR %{y:.2f}<extra></extra>'))
        fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode='lines', hoverinfo='skip', showlegend=False,
                                 line=dict(color=P['muted'], width=1, dash='dot')))
        fig.update_xaxes(title='False-positive rate', range=[0, 1])
        fig.update_yaxes(title='True-positive rate', range=[0, 1.02], scaleanchor='x')
        fig.update_layout(showlegend=False)
        show(styled(fig, f"Holdout ROC · n={h['n_holdout']} ({h['n_pos']} high)", height=400))
        st.dataframe(pd.DataFrame([dict(When=time.strftime('%d %b %H:%M', time.localtime(x['created'])),
                                        Run=x['run_id'][:6], AUC=x['auc'], CI=f"{x['ci_low']:.3f}–{x['ci_high']:.3f}",
                                        Baseline=x['baseline_auc']) for x in peeks]),
                     hide_index=True, width='stretch')

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
with tabs[5]:
    if runs_all:
        flat = pd.DataFrame([{
            'run_id': r['run_id'], 'target': r['target'], 'model': r['model'], 'fs_method': r['fs_method'],
            'auc_mean': r['auc_mean'], 'auc_sd': r['auc_sd'], 'baseline_auc': r['baseline_auc_mean'],
            'gain': r['gain'], 'n_features_mean': r['n_features_mean'],
            'perm_p': perm_summary(r)[0], 'n_perm': perm_summary(r)[1],
            'holdout_auc': (latest_holdout(r['run_id']) or {}).get('auc'),
            'cfg': json.dumps(r['cfg'], sort_keys=True)} for r in runs_all])
        st.download_button('Download runs (CSV)', flat.to_csv(index=False), 'ml_workbench_runs_v2.csv', 'text/csv')
    if focus_run:
        st.download_button('Download focused configuration (JSON)', json.dumps(focus_run['cfg'], indent=2),
                           f"workbench_cfg_{focus_run['run_id']}.json", 'application/json')
    st.caption(f'Run log: {E.RUN_LOG}  ·  fold cache: {E.CACHE_DIR}  ·  Optuna journal: {E.OPTUNA_JOURNAL}')
