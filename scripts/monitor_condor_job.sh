#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 <cluster-id[.proc-id]> [interval-seconds]" >&2
    echo "Example: $0 143159 600" >&2
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage
    exit 2
fi

job_ref=$1
interval=${2:-600}

if [[ ! $job_ref =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Invalid Condor job id: $job_ref" >&2
    exit 2
fi
if [[ ! $interval =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid polling interval: $interval" >&2
    exit 2
fi

status_name() {
    case "$1" in
        1) echo IDLE ;;
        2) echo RUNNING ;;
        3) echo REMOVED ;;
        4) echo COMPLETED ;;
        5) echo HELD ;;
        6) echo TRANSFERRING_OUTPUT ;;
        7) echo SUSPENDED ;;
        *) echo UNKNOWN ;;
    esac
}

timestamp() {
    date '+%Y-%m-%d %H:%M:%S %Z'
}

while :; do
    queue_lines=$(condor_q "$job_ref" -af ClusterId ProcId JobStatus JobStatusText 2>/dev/null || true)
    if [[ -n $queue_lines ]]; then
        printf '[%s] queue status:\n%s\n' "$(timestamp)" "$queue_lines"
        terminal=0
        while read -r cluster proc status status_text; do
            [[ -n ${cluster:-} ]] || continue
            name=$(status_name "${status:-0}")
            printf '[%s] %s.%s %s (%s)\n' "$(timestamp)" "$cluster" "$proc" "$name" "${status_text:-}"
            case "${status:-0}" in
                3|4|5) terminal=1 ;;
            esac
        done <<< "$queue_lines"
        (( terminal == 1 )) && exit 0
    else
        history_lines=$(condor_history "$job_ref" -limit 1 -af ClusterId ProcId JobStatus ExitCode 2>/dev/null || true)
        if [[ -n $history_lines ]]; then
            printf '[%s] history status:\n%s\n' "$(timestamp)" "$history_lines"
            read -r cluster proc status exit_code <<< "$history_lines"
            printf '[%s] %s.%s %s (exit=%s)\n' "$(timestamp)" "$cluster" "$proc" "$(status_name "${status:-0}")" "${exit_code:-unknown}"
            exit 0
        fi
        printf '[%s] %s is not currently visible in queue or history; retrying\n' "$(timestamp)" "$job_ref"
    fi
    sleep "$interval"
done
