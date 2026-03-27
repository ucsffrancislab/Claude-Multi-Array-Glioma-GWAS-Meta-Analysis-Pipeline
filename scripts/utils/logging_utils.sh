#!/usr/bin/env bash
# =============================================================================
# logging_utils.sh — Consistent timestamped logging for all bash scripts
# Source this file: source "$(dirname "$0")/utils/logging_utils.sh"
# =============================================================================

# Colors (disabled if not a terminal)
if [[ -t 1 ]]; then
    _CLR_RED='\033[0;31m'
    _CLR_YLW='\033[0;33m'
    _CLR_GRN='\033[0;32m'
    _CLR_CYN='\033[0;36m'
    _CLR_RST='\033[0m'
    _CLR_BLD='\033[1m'
else
    _CLR_RED='' _CLR_YLW='' _CLR_GRN='' _CLR_CYN='' _CLR_RST='' _CLR_BLD=''
fi

_timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

log_info()    { echo -e "$(_timestamp) ${_CLR_GRN}[INFO]${_CLR_RST}  $*"; }
log_warn()    { echo -e "$(_timestamp) ${_CLR_YLW}[WARN]${_CLR_RST}  $*" >&2; }
log_error()   { echo -e "$(_timestamp) ${_CLR_RED}[ERROR]${_CLR_RST} $*" >&2; }
log_step()    { echo -e "\n$(_timestamp) ${_CLR_CYN}${_CLR_BLD}========== $* ==========${_CLR_RST}"; }
log_substep() { echo -e "$(_timestamp) ${_CLR_CYN}--- $* ---${_CLR_RST}"; }
log_cmd()     { log_info "Running: $*"; "$@"; }

# Print a separator line
log_separator() {
    echo -e "$(_timestamp) ${_CLR_CYN}$(printf '=%.0s' {1..60})${_CLR_RST}"
}

# Log start/end of a timed section
# Usage:  log_timer_start "description"
#         ... do work ...
#         log_timer_end "description"
declare -A _TIMERS 2>/dev/null || true

log_timer_start() {
    local label="$1"
    _TIMERS["$label"]=$(date +%s)
    log_info "Started: ${label}"
}

log_timer_end() {
    local label="$1"
    local start=${_TIMERS["$label"]:-$(date +%s)}
    local elapsed=$(( $(date +%s) - start ))
    local mins=$(( elapsed / 60 ))
    local secs=$(( elapsed % 60 ))
    log_info "Completed: ${label} (${mins}m ${secs}s)"
}

# Die with error message
die() {
    log_error "$@"
    exit 1
}

# Check that a command exists
require_cmd() {
    local cmd="$1"
    command -v "$cmd" &>/dev/null || die "Required command not found: ${cmd}"
}
