#!/usr/bin/env bash
# Idempotent launcher for the local FASH trajectory explorer.
#
#   run_explorer.sh [start|status|stop|restart|log]
#
# start (default) reuses an already-running instance instead of failing on
# "address already in use", waits until the app actually answers, and only then
# opens the browser. The workflowr index page links to the default port, so
# leaving FASH_EXPLORER_PORT unset is what keeps that link working.
#
# Environment:
#   FASH_EXPLORER_PORT   port to serve on (default 7421)
#   FASH_EXPLORER_WAIT   seconds to wait for startup (default 180)

set -uo pipefail

PORT="${FASH_EXPLORER_PORT:-7421}"
WAIT_SECONDS="${FASH_EXPLORER_WAIT:-180}"
URL="http://127.0.0.1:${PORT}/"
LOG_FILE="${TMPDIR:-/tmp}/fash_explorer_${PORT}.log"

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workflowr_root="$(cd "${script_directory}/../.." && pwd)"
if [ ! -f "${workflowr_root}/_workflowr.yml" ]; then
  echo "error: could not find _workflowr.yml above ${script_directory}." >&2
  exit 1
fi
cd "${workflowr_root}" || exit 1

# ---------------------------------------------------------------------------

serving() {
  [ "$(curl -s -m 3 -o /dev/null -w '%{http_code}' "${URL}" 2>/dev/null)" = "200" ]
}

port_holder() {
  lsof -nP -iTCP:"${PORT}" -sTCP:LISTEN -t 2>/dev/null | head -1
}

# Identify the app by who holds the port, then confirm with the command line.
# Do not pattern-match "port = ${PORT}" through pgrep or ps: macOS renders
# spaces inside a single argv element as "~+~", so such a pattern never matches.
app_pid() {
  local pid
  pid="$(port_holder)"
  [ -n "${pid}" ] || return 1
  if ps -o command= -p "${pid}" 2>/dev/null | grep -q "fash_trajectory_explorer"; then
    echo "${pid}"
    return 0
  fi
  return 1
}

open_browser() {
  if command -v open >/dev/null 2>&1; then
    open "${URL}"
  elif command -v xdg-open >/dev/null 2>&1; then
    xdg-open "${URL}"
  else
    echo "Open ${URL} in your browser."
  fi
}

do_status() {
  local holder
  holder="$(port_holder)"
  if serving; then
    echo "running: ${URL} (pid ${holder:-unknown})"
    return 0
  fi
  if [ -n "${holder}" ]; then
    echo "port ${PORT} is held by pid ${holder} but is not answering as the explorer."
    return 2
  fi
  echo "not running: nothing is listening on port ${PORT}."
  return 1
}

do_stop() {
  local pid
  pid="$(app_pid)"
  if [ -z "${pid}" ]; then
    pid="$(port_holder)"
    if [ -n "${pid}" ]; then
      echo "port ${PORT} is held by pid ${pid}, which is not this app. Leaving it alone."
      return 1
    fi
    echo "nothing to stop on port ${PORT}."
    return 0
  fi
  echo "stopping explorer (pid ${pid})..."
  kill "${pid}" 2>/dev/null
  local waited=0
  while kill -0 "${pid}" 2>/dev/null && [ "${waited}" -lt 20 ]; do
    sleep 1
    waited=$((waited + 1))
  done
  if kill -0 "${pid}" 2>/dev/null; then
    echo "pid ${pid} did not exit after ${waited}s; sending SIGKILL."
    kill -9 "${pid}" 2>/dev/null
  fi
  echo "stopped."
}

do_start() {
  if serving; then
    echo "already running on ${URL} (pid $(port_holder)); reusing it."
    open_browser
    return 0
  fi

  local holder
  holder="$(port_holder)"
  if [ -n "${holder}" ]; then
    echo "error: port ${PORT} is held by pid ${holder} but does not answer as the explorer." >&2
    echo "       Stop that process, or set FASH_EXPLORER_PORT to a free port." >&2
    return 1
  fi

  echo "starting the explorer on port ${PORT} (loads the IWP1 fit, expect ~20s)..."
  : > "${LOG_FILE}"
  nohup Rscript --vanilla -e \
    "shiny::runApp('apps/fash_trajectory_explorer', host = '127.0.0.1', port = ${PORT}, launch.browser = FALSE)" \
    >> "${LOG_FILE}" 2>&1 &
  local pid=$!
  # nohup already detaches the child; disown just keeps this shell quiet on exit.
  disown %% 2>/dev/null || true

  local waited=0
  while [ "${waited}" -lt "${WAIT_SECONDS}" ]; do
    if serving; then
      echo "ready after ${waited}s: ${URL}"
      echo "log: ${LOG_FILE}"
      open_browser
      return 0
    fi
    if ! kill -0 "${pid}" 2>/dev/null; then
      echo "error: the R process exited during startup. Last lines of ${LOG_FILE}:" >&2
      tail -20 "${LOG_FILE}" >&2
      return 1
    fi
    sleep 1
    waited=$((waited + 1))
  done

  echo "error: the app did not answer within ${WAIT_SECONDS}s. Last lines of ${LOG_FILE}:" >&2
  tail -20 "${LOG_FILE}" >&2
  return 1
}

case "${1:-start}" in
  start) do_start ;;
  status) do_status ;;
  stop) do_stop ;;
  restart)
    do_stop
    do_start
    ;;
  log)
    if [ -f "${LOG_FILE}" ]; then
      tail -40 "${LOG_FILE}"
    else
      echo "no log at ${LOG_FILE} yet."
    fi
    ;;
  *)
    echo "usage: $(basename "$0") [start|status|stop|restart|log]" >&2
    exit 1
    ;;
esac
