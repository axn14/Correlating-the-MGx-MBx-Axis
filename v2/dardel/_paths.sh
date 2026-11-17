# Shared helpers for the laptop-side scripts. Sourced, not executed.
#
# Two problems, both specific to driving Dardel from Windows.
#
# 1. The Windows-side shells mount drives at different roots, and dardel.env has to
#    work in whichever one you happen to be in.
# 2. The PDC key is usually not at a default identity name (MobaXterm names it
#    id-ed25519-pdc), so ssh will not pick it up on its own. configure_ssh makes
#    every ssh and rsync call in these scripts use PDC_SSH_KEY when it is set.
#
# The drive roots:
#
#   WSL          /mnt/d/KI Research
#   MobaXterm    /drives/d/KI Research
#   Git Bash     /d/KI Research
#
# resolve_local_path takes a path in any of those three forms and returns the one that
# exists here, so the same dardel.env works from all three without editing.

resolve_local_path() {
    local p="$1" letter rest cand
    [ -d "$p" ] && { printf '%s' "$p"; return 0; }
    case "$p" in
        /mnt/?/*)    letter="${p:5:1}"; rest="${p:7}" ;;
        /drives/?/*) letter="${p:8:1}"; rest="${p:10}" ;;
        /?/*)        letter="${p:1:1}"; rest="${p:3}" ;;
        *)           printf '%s' "$p"; return 0 ;;
    esac
    for cand in "/mnt/${letter}/${rest}" "/drives/${letter}/${rest}" "/${letter}/${rest}"; do
        [ -d "$cand" ] && { printf '%s' "$cand"; return 0; }
    done
    printf '%s' "$p"
}

# LOCAL_PULL usually does not exist yet, so it is translated by prefix rather than by
# testing for the directory: whichever drive root LOCAL_ROOT resolved to wins.
retarget_to_same_root() {
    local p="$1" root_resolved="$2" root_given="$3"
    local prefix_resolved prefix_given
    case "$root_resolved" in
        /mnt/?/*)    prefix_resolved="/mnt/${root_resolved:5:1}" ;;
        /drives/?/*) prefix_resolved="/drives/${root_resolved:8:1}" ;;
        /?/*)        prefix_resolved="/${root_resolved:1:1}" ;;
        *)           printf '%s' "$p"; return 0 ;;
    esac
    case "$root_given" in
        /mnt/?/*)    prefix_given="/mnt/${root_given:5:1}" ;;
        /drives/?/*) prefix_given="/drives/${root_given:8:1}" ;;
        /?/*)        prefix_given="/${root_given:1:1}" ;;
        *)           printf '%s' "$p"; return 0 ;;
    esac
    [ "$prefix_resolved" = "$prefix_given" ] && { printf '%s' "$p"; return 0; }
    case "$p" in
        "${prefix_given}"/*) printf '%s%s' "$prefix_resolved" "${p#"$prefix_given"}" ;;
        *) printf '%s' "$p" ;;
    esac
}

# Call after sourcing dardel.env. Rewrites LOCAL_ROOT and LOCAL_PULL for this shell.
normalise_local_paths() {
    local given="${LOCAL_ROOT}"
    LOCAL_ROOT="$(resolve_local_path "$LOCAL_ROOT")"
    if [ -n "${LOCAL_PULL:-}" ]; then
        LOCAL_PULL="$(retarget_to_same_root "$LOCAL_PULL" "$LOCAL_ROOT" "$given")"
    fi
    if [ ! -d "$LOCAL_ROOT" ]; then
        echo "LOCAL_ROOT not found in this shell: ${given}" >&2
        echo "Tried /mnt, /drives and /<letter> drive roots. Shells seen here:" >&2
        for r in /mnt /drives; do [ -d "$r" ] && echo "  ${r}: $(ls -1 "$r" 2>/dev/null | tr '\n' ' ')" >&2; done
        return 1
    fi
    return 0
}

# Route every ssh and rsync call through PDC_SSH_KEY when dardel.env sets one.
# The ssh function shadows the binary for direct calls; RSYNC_RSH does the same for
# rsync, which spawns its own remote shell and would not see a shell function.
# Nothing here reads or copies the key: ssh does, as it normally would.
configure_ssh() {
    [ -n "${PDC_SSH_KEY:-}" ] || return 0
    local key
    key="$(resolve_local_path "$PDC_SSH_KEY")"
    # resolve_local_path tests for a directory, so for a file fall back to the literal
    # value and let ssh report it if it is wrong.
    [ -f "$key" ] || key="$PDC_SSH_KEY"
    if [ ! -f "$key" ]; then
        echo "PDC_SSH_KEY is set but no such file: ${PDC_SSH_KEY}" >&2
        echo "Leave it empty to use ssh's own default identities or agent." >&2
        return 1
    fi
    # A global, not the local: the function body is evaluated when it is CALLED, by
    # which point a local would have gone out of scope and -i would get an empty path.
    PDC_SSH_KEY_RESOLVED="$key"
    ssh() { command ssh -i "$PDC_SSH_KEY_RESOLVED" "$@"; }
    export RSYNC_RSH="ssh -i '${key}'"
    echo "== ssh identity: ${key}"
    return 0
}
