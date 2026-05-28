from __future__ import annotations
import streamlit as st
import pandas as pd
from PIL import Image
from config import FIGURES_DIR, TABLES_DIR


def load_figure(stem: str) -> Image.Image | None:
    """Return the first image in FIGURES_DIR whose name starts with stem."""
    matches = sorted(FIGURES_DIR.glob(f"{stem}*"))
    if not matches:
        st.warning(f"Figure not found: {stem}* in {FIGURES_DIR}")
        return None
    return Image.open(matches[0])


def show_figure(stem: str, caption: str = "", use_container_width: bool = True) -> None:
    img = load_figure(stem)
    if img is not None:
        st.image(img, caption=caption, use_container_width=use_container_width)


def load_csv(name: str) -> pd.DataFrame | None:
    """Load a CSV — prefers user-run results if available in session state."""
    import streamlit as st
    user_dir = getattr(st.session_state, "run_output_dir", None)
    if user_dir is not None:
        user_path = Path(user_dir) / "tables" / name
        if user_path.exists():
            try:
                return pd.read_csv(user_path)
            except Exception:
                pass
    path = TABLES_DIR / name
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def metric_row(items: list[tuple[str, str]]) -> None:
    """Render (label, value) pairs as a horizontal row of st.metric cards."""
    cols = st.columns(len(items))
    for col, (label, value) in zip(cols, items):
        col.metric(label, value)
