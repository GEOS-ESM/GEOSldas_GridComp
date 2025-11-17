#!/usr/bin/env bash
set -euo pipefail

die(){ echo "ERROR: $*" >&2; exit 1; }

# --- self-locate experiment root & ids; allow env/.env overrides ---
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../<EXPID>/regress
expdir_default="$(cd "$script_dir/.." && pwd)"               # -> <EXPDIR>
expid_default="$(basename "$expdir_default")"

# optional overrides from a local .env in regress/
[[ -f "$script_dir/.env" ]] && . "$script_dir/.env"

# --- locate regression home and templates ---
REGRESS_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE_DIR="${REGRESS_TEMPLATE_DIR:-$REGRESS_HOME/templates}"

# if installed under share/geosldas/regression (optional)
if [[ ! -d "$TEMPLATE_DIR" && -n "${GEOSLDAS_SHARE:-}" ]]; then
  TEMPLATE_DIR="$GEOSLDAS_SHARE/regression/templates"
fi

EXPDIR="${EXPDIR:-$expdir_default}"
EXPID="${EXPID:-$expid_default}"

# Try to infer EXPDOMAIN from existing output/<DOMAIN>/ if not provided
if [[ -z "${EXPDOMAIN:-}" && -d "$EXPDIR/output" ]]; then
  EXPDOMAIN="$(find "$EXPDIR/output" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | head -n1 || true)"
fi

[[ -d "$EXPDIR/run" ]] || die "Cannot find run/ under EXPDIR=$EXPDIR"
: "${EXPDOMAIN:?set EXPDOMAIN=your_domain (or create $EXPDIR/output/<DOMAIN>)}"

SUBMIT="${SUBMIT:-sbatch}"
[[ "$SUBMIT" == "sbatch" ]] || die "This regression driver supports Slurm only (SUBMIT=sbatch)."
JOBFILE="${JOBFILE:-$EXPDIR/run/lenkf.j}"

# --- Optional layout-invariance test (off by default) ---
RUN_LAYOUT="${RUN_LAYOUT:-0}"          # enable with RUN_LAYOUT=1
ABS_TOL="${ABS_TOL:-1e-15}"            # nccmp absolute tolerance (-t)
REL_TOL="${REL_TOL:-1e-12}"             # nccmp relative percent tolerance (-T)
NCCMP_FLAGS_TOL="${NCCMP_FLAGS_TOL:--dmfMNS -G history -t $ABS_TOL -T $REL_TOL}"

# --- ensure compare tools available ---
module load nccmp || echo "(warn) nccmp not found in module path"

NCCMP="${NCCMP:-nccmp}"

T1_HHMMSS="240000"
T2_HHMMSS="120000"
T3_HHMMSS="120000"

# Always patch sandbox HISTORY.rc to 6-hour averages (060000) for regression
HIST_MODE="6h"
# Step size in seconds for HISTORY collection (default 6h = 21600)
HIST_STEP_SEC="${HIST_STEP_SEC:-21600}"
HIST_STEP_OFFSET_SEC="${HIST_STEP_OFFSET_SEC:-10800}"  # +3 hours (centered timestamps: 03/09/15/21Z)

ROOT="$EXPDIR/regress"
LOGS="$ROOT/logs"
SETS="$ROOT/sets"
mkdir -p "$LOGS" "$SETS"

# --- save on screen log to dir for records
ts=$(date +%Y%m%d_%H%M%S)
exec > >(tee -a "$LOGS/run_${ts}.log") 2>&1
echo "Run started $(date)"

# --- sandbox mode (isolate regression from the real experiment) ---
SANDBOX_ROOT="$EXPDIR/regress/sandbox"
SANDBOX_EXPDIR="$SANDBOX_ROOT/$EXPID"     # e.g., .../regress/sandbox/CURRENT
SANDBOX_OUT="$SANDBOX_EXPDIR/output/$EXPDOMAIN"

make_sandbox() {
  echo "Creating regression sandbox at: $SANDBOX_EXPDIR"
  mkdir -p "$SANDBOX_EXPDIR"

  #  Copy run/ (small files)
  rsync -a --delete "$EXPDIR/run/"   "$SANDBOX_EXPDIR/run/"

  # Lightweight input/: symlink everything from real input/* except restart
  mkdir -p "$SANDBOX_EXPDIR/input"
  for e in "$EXPDIR/input/"* ; do
    bn="$(basename "$e")"
    [[ "$bn" == "restart" ]] && continue
    ln -sfn "$e" "$SANDBOX_EXPDIR/input/$bn"
  done

  #  Build link to build/ 
  ln -sfn "$EXPDIR/build" "$SANDBOX_EXPDIR/build"

  # Create scratch/ and sandbox OUT path
  mkdir -p "$SANDBOX_EXPDIR/scratch" "$SANDBOX_OUT"

  # Patch OUT_PATH in sandbox LDAS.rc to sandbox output
  sed -i -E "s|^OUT_PATH:\s*.*$|OUT_PATH:                         $SANDBOX_OUT|" \
      "$SANDBOX_EXPDIR/run/LDAS.rc"

  # Patch EXPDIR in sandbox lenkf.j to sandbox path
  sed -i -E "s|^(setenv[[:space:]]+EXPDIR[[:space:]]+).*|\1$SANDBOX_EXPDIR|" "$SANDBOX_EXPDIR/run/lenkf.j"

  # Force the job to use the sandbox HISTORY.rc and enable history 
  {
    echo "setenv HISTORY_RC $SANDBOX_EXPDIR/run/HISTORY.rc"
    echo "setenv DO_HISTORY TRUE"
    echo "setenv DO_HIST TRUE"
  } >> "$SANDBOX_EXPDIR/run/lenkf.j"

  # keep POSTPROC_HIST=0 so 6-hour sub-daily files are not bundled away
  sed -i -E 's/^(setenv[[:space:]]+POSTPROC_HIST[[:space:]]+).*/\10/' "$SANDBOX_EXPDIR/run/lenkf.j"

  # also patch any wrapper hard codes
  if [[ -f "$SANDBOX_EXPDIR/run/globalcs.model.exec" ]]; then
    sed -i -E "s|^(HISTRC_FILE:\s*).*$|\\1$SANDBOX_EXPDIR/run/HISTORY.rc|" "$SANDBOX_EXPDIR/run/globalcs.model.exec"
  fi
  sed -i -E "s|^( *HIST_CF:\s*).*$|\\1$SANDBOX_EXPDIR/run/HISTORY.rc|" "$SANDBOX_EXPDIR/run/CAP.rc"

  # Patch job log destinations to sandbox scratch (slurm logs)
  sed -i -E "s|^#SBATCH --output=.*|#SBATCH --output=$SANDBOX_EXPDIR/scratch/GEOSldas_log_txt|" \
      "$SANDBOX_EXPDIR/run/lenkf.j"
  sed -i -E "s|^#SBATCH --error=.*|#SBATCH --error=$SANDBOX_EXPDIR/scratch/GEOSldas_err_txt|" \
      "$SANDBOX_EXPDIR/run/lenkf.j"

  # Copy ONLY the original start restarts into sandbox rs/
  local yy="${ORIG_NYMD:0:4}" mm="${ORIG_NYMD:4:2}"
  local srs="$SANDBOX_OUT/rs/ens0000/Y$yy/M$mm"
  mkdir -p "$srs"

  for comp in catch landice; do
    src="$EXPDIR/output/$EXPDOMAIN/rs/ens0000/Y$yy/M$mm/$EXPID.${comp}_internal_rst.${ORIG_NYMD}_${ORIG_NHMS}"
    [[ -f "$src" ]] || src="$EXPDIR/output/$EXPDOMAIN/rs/ens0000/Y$yy/M$mm/$EXPID.${comp}_internal_rst.${ORIG_NYMD}_${ORIG_NHMS:0:4}"
    if [[ -f "$src" ]]; then
      cp -p "$src" "$srs/"
    else
      echo "(warn) missing original $comp restart at ${ORIG_NYMD}_${ORIG_NHMS} (HHMMSS/HHMM) in real rs/"
    fi
  done

  #  Make sandbox input/restart/ that points to the **sandbox** copies
  mkdir -p "$SANDBOX_EXPDIR/input/restart"
  for comp in catch landice; do
    tgt="$(ls "$srs/$EXPID.${comp}_internal_rst.${ORIG_NYMD}"_* 2>/dev/null | head -n1 || true)"
    if [[ -n "$tgt" ]]; then
      ln -rsf "$tgt" "$SANDBOX_EXPDIR/input/restart/${comp}_internal_rst"
    fi
  done
  # Also bring vegdyn restart local so all three land restarts are inside sandbox
  local vsrc="$EXPDIR/output/$EXPDOMAIN/rs/ens0000/$EXPID.vegdyn_internal_rst"
  if [[ -f "$vsrc" ]]; then
    mkdir -p "$SANDBOX_OUT/rs/ens0000"
    cp -p "$vsrc" "$SANDBOX_OUT/rs/ens0000/"
    ln -sfn "$SANDBOX_OUT/rs/ens0000/$(basename "$vsrc")" \
            "$SANDBOX_EXPDIR/input/restart/vegdyn_internal_rst"
    echo "Linked vegdyn_internal_rst -> $SANDBOX_EXPDIR/input/restart/vegdyn_internal_rst"
  else
    echo "(warn) vegdyn_internal_rst not found at $vsrc; LDAS.rc may still require it"
  fi
  # Bring over the full rc_out so the model can read/write locally
  local rc_real="$EXPDIR/output/$EXPDOMAIN/rc_out"
  local rc_sb="$SANDBOX_OUT/rc_out"
  mkdir -p "$SANDBOX_OUT"
  echo "Syncing rc_out from $rc_real -> $rc_sb ..."
  rsync -a --delete "$rc_real/" "$rc_sb/"

  # Ensure the month directories we’re about to write into exist
  mkdir -p "$rc_sb/Y$yy/M$mm"
  local ny="$yy" nm=$((10#$mm + 1))
  if (( nm > 12 )); then nm=1; ny=$((10#$yy + 1)); fi
  printf -v nm "%02d" "$nm"
  mkdir -p "$rc_sb/Y$ny/M$nm"

  # Now that rc_out is synced, link mwRTM param if present
  shopt -s nullglob
  for m in "$SANDBOX_OUT/rc_out"/Y*/M*/"$EXPID".ldas_mwRTMparam.*.nc4; do
    ln -sfn "$m" "$SANDBOX_EXPDIR/input/restart/mwrtm_param_rst"
    echo "Linked mwrtm_param_rst -> $m"
    break
  done
  shopt -u nullglob
  # Ensure LDAS.rc uses sandbox HISTORY.rc and HISTORY is enabled (if keys exist)
  if grep -qE '^\s*HISTORY(_RC)?\s*:' "$SANDBOX_EXPDIR/run/LDAS.rc"; then
    sed -i -E 's|^\s*HISTORY(_RC)?\s*:\s*.*$|HISTORY_RC:                     HISTORY.rc|' "$SANDBOX_EXPDIR/run/LDAS.rc"
  fi
  if grep -qE '^\s*(DO_HIST|DO_HISTORY)\s*:' "$SANDBOX_EXPDIR/run/LDAS.rc"; then
    sed -i -E 's|^\s*(DO_HIST|DO_HISTORY)\s*:.*$|\1:                        .true.|' "$SANDBOX_EXPDIR/run/LDAS.rc"
  fi

  # Detect whether landice restart exists in sandbox
  if [[ -L "$SANDBOX_EXPDIR/input/restart/landice_internal_rst" || \
        -f "$SANDBOX_EXPDIR/input/restart/landice_internal_rst" ]]; then
    HAS_LANDICE=1
    echo "Sandbox: landice_internal_rst detected (HAS_LANDICE=1)"
  else
    HAS_LANDICE=0
    echo "Sandbox: NO landice_internal_rst (HAS_LANDICE=0)"
  fi

  echo "Sandbox ready."
}

# -------- helpers -------------------------------------------------------------
patch_history_6h() {
  local hrc="$EXPDIR/run/HISTORY.rc"
  [[ -f "$hrc" ]] || { echo "(warn) HISTORY.rc not found at $hrc"; return 0; }
  echo "Patching HISTORY.rc for 6-hour tavg (060000)…"

  # CF (2d/Nx)
  sed -i -E \
    -e 's|^([[:space:]]*tavg06_2d_lfs_Nx\.frequency:[[:space:]]*)[0-9]{6}|\1060000|' \
    -e 's|^([[:space:]]*tavg06_2d_lnd_Nx\.frequency:[[:space:]]*)[0-9]{6}|\1060000|' \
    -e 's|^([[:space:]]*tavg06_2d_glc_Nx\.frequency:[[:space:]]*)[0-9]{6}|\1060000|' \
    "$hrc"

  # EASE (1d/Nt)
  sed -i -E \
    -e 's|^([[:space:]]*tavg06_1d_lfs_Nt\.frequency:[[:space:]]*)[0-9]{6}|\1060000|' \
    -e 's|^([[:space:]]*tavg06_1d_lnd_Nt\.frequency:[[:space:]]*)[0-9]{6}|\1060000|' \
    -e 's|^([[:space:]]*tavg06_1d_glc_Nt\.frequency:[[:space:]]*)[0-9]{6}|\1060000|' \
    "$hrc"

  # Show what we have now (both Nx and Nt if present)
  grep -E 'tavg06_(1d_.*_Nt|2d_.*_Nx)\.frequency:' "$hrc" | sed 's/^/  /' || true
}

# Write cap_restart to a given NYMD NHMS
write_cap_restart() {
  local ny="$1" nh="$2"
  printf "%s %s\n" "$ny" "$nh" > "$EXPDIR/run/cap_restart"
}

# Build a path to a component restart at a given timestamp (ens0000; adjust if needed)
restart_path_for_stamp() {
  local comp="$1" ny="$2" nh="$3"
  local yy=${ny:0:4} mm=${ny:4:2}
  local base="$EXPDIR/output/$EXPDOMAIN/rs/ens0000/Y${yy}/M${mm}/${EXPID}.${comp}_internal_rst"
  # prefer HHMMSS
  if [[ -f "${base}.${ny}_${nh}" ]]; then
    printf "%s.%s_%s" "$base" "$ny" "$nh"
  elif [[ -f "${base}.${ny}_${nh:0:4}" ]]; then
    printf "%s.%s_%s" "$base" "$ny" "${nh:0:4}"
  else
    # nothing; return empty (caller will warn)
    printf ""
  fi
}

# Repoint input/restart symlinks and reset cap_restart back to the ORIGINAL start
restore_initial_restarts() {
  local ny="$ORIG_NYMD" nh="$ORIG_NHMS"
  echo "Restoring input/restart links and cap_restart to ${ny} ${nh} ..."
  write_cap_restart "$ny" "$nh"

  # Catch and landice have timestamped files; vegdyn is timeless here
  for comp in catch landice; do
    local tgt; tgt="$(restart_path_for_stamp "$comp" "$ny" "$nh")"
    if [[ -n "$tgt" && -f "$tgt" ]]; then
      ln -rsf "$tgt" "$EXPDIR/input/restart/${comp}_internal_rst"
      echo "  -> ${comp}_internal_rst -> $(readlink -f "$EXPDIR/input/restart/${comp}_internal_rst")"
    else
      echo "  (warn) missing expected $comp restart at ${ny}_${nh} (tried HHMMSS and HHMM)"
    fi
  done  # 

  # vegdyn: keep existing link if present (no timestamp)
  if [[ -L "$EXPDIR/input/restart/vegdyn_internal_rst" || -f "$EXPDIR/input/restart/vegdyn_internal_rst" ]]; then
    echo "  -> vegdyn_internal_rst unchanged ($(readlink -f "$EXPDIR/input/restart/vegdyn_internal_rst" 2>/dev/null || true))"
  fi
}

cap_restart_start() {
  # prints: NYMD NHMS
  local f="$EXPDIR/run/cap_restart"
  [[ -s "$f" ]] || die "Missing $f"
  awk '{print $1,$2}' "$f" | head -1
}

set_cap_rc() {
  # Args: JOB_SGMT_HHMMSS END_NYMD END_NHMS
  local job="$1" endy="$2" endh="$3" cap="$EXPDIR/run/CAP.rc"
  [[ -f "$cap" ]] || die "Missing $cap"
  sed -i -E "s/^JOB_SGMT:\s+[0-9]{8}\s+[0-9]{6}/JOB_SGMT:  00000000 ${job}/" "$cap"
  sed -i -E "s/^END_DATE:\s+[0-9]{8}\s+[0-9]{6}/END_DATE:  ${endy} ${endh}/" "$cap"
}

submit_and_wait_done() {
  local out jid state
  out="$($SUBMIT "$JOBFILE")" || die "submit failed"
  echo "Submitted: $out"

  # parse 'Submitted batch job 12345'
  jid="$(echo "$out" | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) {print $i; exit}}')"
  echo "$jid" > "$EXPDIR/scratch/LAST_JOB_ID" || true
  [[ -n "$jid" ]] || die "could not parse Slurm job id from: $out"

  # wait until job leaves the queue
  local tries=800
  while ((tries--)); do
    if [[ -n "$(squeue -h -j "$jid" 2>/dev/null)" ]]; then
      sleep 10; continue
    fi
    sleep 5; break
  done

  # check final state and print sandbox logs if failed
  if command -v sacct >/dev/null 2>&1; then
    state="$(sacct -j "$jid" --format=State%20,ExitCode%10 -P -n 2>/dev/null | head -n1)"
    echo "Job $jid state: $state"
    if echo "$state" | grep -qE 'FAILED|CANCELLED|TIMEOUT'; then
      errf="$EXPDIR/scratch/GEOSldas_err_txt"
      logf="$EXPDIR/scratch/GEOSldas_log_txt"
      echo "---- sandbox error (head) $errf ----"; [[ -f "$errf" ]] && sed -n '1,120p' "$errf" || echo "(no err file)"
      echo "---- sandbox log   (head) $logf ----"; [[ -f "$logf" ]] && sed -n '1,60p'  "$logf" || echo "(no log file)"
      die "batch job $jid failed: $state"
    fi
  fi
}

wait_for_restart_stamp() {
  local y h dir h4 tries
  y="$1"
  h="$2"
  dir="$EXPDIR/output/$EXPDOMAIN/rs"
  h4="${h:0:4}"
  tries=180
  echo "Waiting for restarts at ${y}_${h} (or ${y}_${h4}) in $dir ..."
  while ((tries--)); do
    if find "$dir" -type f \( -name "${EXPID}.*.${y}_${h}*" -o -name "${EXPID}.*.${y}_${h4}*" \) | grep -q .; then
      echo "Found restarts:"
      find "$dir" -type f \( -name "${EXPID}.*.${y}_${h}*" -o -name "${EXPID}.*.${y}_${h4}*" \) -ls
      return 0
    fi
    sleep 5
  done
  die "No restarts found for ${y}_${h}/${y}_${h4} in $dir"
}
collect_set() {
  local tag y h tgt h4 dir
  tag="$1"
  y="$2"
  h="$3"
  tgt="$SETS/$tag"
  h4="${h:0:4}"
  dir="$EXPDIR/output/$EXPDOMAIN/rs"
  mkdir -p "$tgt"
  echo "Collecting into $tgt ..."
  # List both HHMMSS and HHMM matches, then copy each
  while IFS= read -r f; do
    [[ -n "$f" ]] || continue
    echo "  -> $f"
    cp -p "$f" "$tgt/"
  done < <(find "$dir" -type f \( -name "${EXPID}.*.${y}_${h}*" -o -name "${EXPID}.*.${y}_${h4}*" \) | sort)
}

compare_sets() {
  # Args: tagA tagB
  local A="$SETS/$1" B="$SETS/$2"
  [[ -d "$A" && -d "$B" ]] || die "Missing compare dirs: $A or $B"

  # Enforce identical file lists (restart files)
  local listA listB
  listA="$(find "$A" -maxdepth 1 -type f -printf '%f\0' | sort -z | tr '\0' '\n')"
  listB="$(find "$B" -maxdepth 1 -type f -printf '%f\0' | sort -z | tr '\0' '\n')"
  if ! diff -u <(printf '%s\n' "$listA") <(printf '%s\n' "$listB") >/dev/null; then
    echo "File list mismatch between $A and $B:"
    diff -u <(printf '%s\n' "$listA") <(printf '%s\n' "$listB") || true
    return 2
  fi

  # One-to-one content compare
  local rc=0
  while IFS= read -r -d '' fa; do
    local base fb
    base="$(basename "$fa")"
    fb="$B/$base"
    if [[ -f "$fb" ]]; then
      nc_compare restart "$fa" "$fb" "$base" || rc=1
    else
      echo "Missing in $B: $base"; rc=1
    fi
  done < <(find "$A" -maxdepth 1 -type f -print0)

  return $rc
}

promote_restarts_from_stamp() {
  local y h h4 link target
  y="$1"
  h="$2"
  h4="${h:0:4}"
  for link in "$EXPDIR/input/restart/"*; do
    [[ -L "$link" || -e "$link" ]] || continue
    case "$(basename "$link")" in
      *internal_rst|*rst)
        target="$(find "$EXPDIR/output/$EXPDOMAIN/rs" -type f \
                  \( -name "${EXPID}.*.${y}_${h}*" -o -name "${EXPID}.*.${y}_${h4}*" \) \
                  | grep -E "$(basename "$link" | sed 's/_internal_rst$//;s/_rst$//')" \
                  | head -n1 || true)"
        if [[ -n "$target" ]]; then
          rm -f "$link"
          ln -rs "$target" "$link"
          echo "Promoted $(basename "$link") -> $(readlink -f "$link")"
        fi
        ;;
    esac
  done
}
# Test History output as well lnd lfs!!
collect_hist_set() {
  # Args: tag yyyymmdd hhmm[ss]
  local tag="$1" y="${2:?collect_hist_set needs yyyymmdd}" h="${3:-}"
  if [[ -z "${h}" ]]; then
    echo "(warn) collect_hist_set: missing HHMM/HHMMSS for ${y}; skipping HISTORY collect for tag=${tag}"
    return 0
  fi

  local h4="${h:0:4}"
  local tgt="$SETS/$tag"; mkdir -p "$tgt"
  local catdir="$EXPDIR/output/$EXPDOMAIN/cat/ens0000"
  echo "Collecting HISTORY into $tgt for ${y}_${h}..."

  # helper: list files for a specific date+stamp
  _list_hist_for_stamp() {
    local yy="$1" hh="$2"
    find "$catdir" -type f -name "${EXPID}.*.${yy}_${hh}z*.nc4" -print0
  }

  # Try exact match first (HHMMSS or HHMM)
  mapfile -d '' files < <(find "$catdir" -type f \
    \( -name "${EXPID}.*.${y}_${h}z*.nc4" \
       -o -name "${EXPID}.*.${y}_${h4}z*.nc4" \) \
    -print0)

  # If none and we asked for 00Z, try previous day 12Z (tavg06 time averaged)
  if (( ${#files[@]} == 0 )) && [[ "${h4}" == 0000 ]]; then
    local TICK="$EXPDIR/build/bin/tick"
    if [[ -x "$TICK" ]]; then
      read prevY prevH < <("$TICK" "$y" "${h4}00" -43200)   # minus 12h
      local prevH4="${prevH:0:4}"
      if [[ "$prevH4" == "1200" ]]; then
        echo "(note) requested ${y}_${h4}, using previous day ${prevY}_1200 (tavg06 cadence)"
        mapfile -d '' files < <(_list_hist_for_stamp "$prevY" "1200")
      fi
    fi
  fi

  # If still none, pick best available stamp for that date:
  #    prefer 12Z tavg, else 00Z inst, else latest
  if (( ${#files[@]} == 0 )); then
    local found_any=0 found_tavg_1200=0 found_inst_0000=0 latest=""
    while IFS= read -r -d '' f; do
      found_any=1
      local b stream ts
      b="$(basename "$f")"
      stream="$(awk -F. '{print $2}' <<< "$b")"
      ts="$(awk -F'[._]' -v d="$y" '{
               for (i=1;i<=NF;i++) if ($i ~ d) {print $(i+1); exit}
             }' <<< "$b")"
      ts="${ts//[^0-9]/}"
      [[ -n "$ts" ]] || continue
      if [[ "$stream" == tavg* && "$ts" == 1200* ]]; then found_tavg_1200=1; fi
      if [[ "$stream" == inst* && "$ts" == 0000* ]]; then found_inst_0000=1; fi
      if [[ -z "$latest" || "$ts" > "$latest" ]]; then latest="$ts"; fi
    done < <(find "$catdir" -type f -name "${EXPID}.*.${y}_*z*.nc4" -print0)

    if (( ! found_any )); then
      echo "(warn) no HISTORY files found for ${y}_${h} (or ${y}_${h4}) in $catdir"
      return 0
    fi

    local chosen="" reason=""
    if (( found_tavg_1200 )); then
      chosen="1200"; reason="tavg06 streams write at 12Z"
    elif (( found_inst_0000 )); then
      chosen="0000"; reason="instantaneous streams write at 00Z"
    else
      chosen="$latest"; reason="using latest available stamp for the date"
    fi
    echo "(note) requested ${y}_${h4}, using available stamp ${y}_${chosen} instead (${reason})"
    mapfile -d '' files < <(_list_hist_for_stamp "$y" "$chosen")
  fi

  # If still none, nothing to do
  if (( ${#files[@]} == 0 )); then
    echo "(warn) no HISTORY files to collect after cadence fallback"
    return 0
  fi

  # Summarize streams
  echo "HISTORY streams in this set:"
  printf "%s\0" "${files[@]}" \
    | xargs -0 -n1 basename \
    | awk -F. '{print $2}' \
    | sort | uniq -c \
    | awk '{printf "  %s=%s\n", $2, $1}'

  # Copy files
  for f in "${files[@]}"; do
    echo "  -> $f"
    cp -p "$f" "$tgt/"
  done
}
compare_hist_sets() {
  # Args: tagA tagB
  local A="$SETS/$1" B="$SETS/$2"
  [[ -d "$A" && -d "$B" ]] || { echo "No history sets to compare ($A,$B)"; return 0; }

  echo "Comparing HISTORY sets: $1 vs $2"

  summarize_hist_dir() {
    local d="$1"
    echo "  Streams in $(basename "$d"):"
    find "$d" -maxdepth 1 -type f -name "*.nc4" -print0 \
      | xargs -0 -n1 basename 2>/dev/null \
      | awk -F. 'NF>=3 {print $2}' \
      | sort | uniq -c \
      | awk '{printf "    %s=%s\n", $2, $1}'
  }
  summarize_hist_dir "$A"
  summarize_hist_dir "$B"

  # If both empty, warn and mark as not comparable
  local cntA cntB
  cntA=$(find "$A" -maxdepth 1 -type f -name "*.nc4" | wc -l)
  cntB=$(find "$B" -maxdepth 1 -type f -name "*.nc4" | wc -l)
  if [[ "$cntA" -eq 0 && "$cntB" -eq 0 ]]; then
    echo "(warn) No HISTORY files in either set; skipping comparison."
    return 1
  fi

  # Enforce identical *.nc4 file lists
  local listA listB
  listA="$(find "$A" -maxdepth 1 -type f -name '*.nc4' -printf '%f\0' | sort -z | tr '\0' '\n')"
  listB="$(find "$B" -maxdepth 1 -type f -name '*.nc4' -printf '%f\0' | sort -z | tr '\0' '\n')"

  if ! diff -u <(printf '%s\n' "$listA") <(printf '%s\n' "$listB") >/dev/null; then
    echo "HISTORY file list mismatch between $A and $B:"
    diff -u <(printf '%s\n' "$listA") <(printf '%s\n' "$listB") || true
    return 2
  fi

  # One-to-one strict compare of *.nc4
  local rc=0
  while IFS= read -r -d '' fa; do
    local base fb
    base="$(basename "$fa")"
    fb="$B/$base"
    if [[ -f "$fb" ]]; then
       nc_compare history "$fa" "$fb" "$base" || rc=1
    else
      echo "HIST MISSING in $(basename "$B"): $base"; rc=1
    fi
  done < <(find "$A" -maxdepth 1 -type f -name '*.nc4' -print0)

  return $rc
}
# Collect hourly HISTORY for (start, end]: centers at +3h, +9h, +15h, +21h, etc.
collect_hist_range() {
  # Args: tag start_yyyymmdd start_hhmmss end_yyyymmdd end_hhmmss
  local tag="$1" sy="$2" sh="$3" ey="$4" eh="$5"
  local tgt="$SETS/$tag"; mkdir -p "$tgt"
  local catdir="$EXPDIR/output/$EXPDOMAIN/cat/ens0000"
  local TICK="$EXPDIR/build/bin/tick"
  [[ -x "$TICK" ]] || { echo "(warn) tick not found; cannot collect hourly HISTORY"; return 0; }

  echo "Collecting 6H HISTORY (centered) into $tgt for ($sy $sh, $ey $eh] ..."

  local cy ch
  read cy ch < <("$TICK" "$sy" "$sh" "$HIST_STEP_OFFSET_SEC")  # first center (e.g., +3h)
  while :; do
    if [[ "$cy$ch" > "$ey$eh" ]]; then break; fi
    local h4="${ch:0:4}"

mapfile -d '' files < <(
  find "$catdir" -type f \
    \( -name "${EXPID}.*.${cy}_${ch}z*.nc4" -o -name "${EXPID}.*.${cy}_${h4}z*.nc4" \) \
  -print0 | sort -z
)
    if (( ${#files[@]} == 0 )); then
      echo "(warn) missing 6H center HISTORY for stamp ${cy}_${ch} (or ${cy}_${h4})"
    else
      # Summarize streams and copy
      printf "%s\0" "${files[@]}" \
        | xargs -0 -n1 basename \
        | awk -F. '{print $2}' \
        | sort | uniq -c \
        | awk -v ts="${cy}_${h4}" '{printf "  [%s] %s=%s\n", ts, $2, $1}'
      for f in "${files[@]}"; do cp -p -- "$f" "$tgt/"; done
    fi

    read cy ch < <("$TICK" "$cy" "$ch" "$HIST_STEP_SEC")       # advance by 6h
  done
}

ensure_catparam_next00z() {
  local rc_sb="$EXPDIR/output/$EXPDOMAIN/rc_out"
  local base="$EXPID.ldas_catparam"
  local sy="$1" sh="$2"
  local TICK="$EXPDIR/build/bin/tick"
  [[ -x "$TICK" ]] || return 0

  local sY=${sy:0:4} sM=${sy:4:2}
  local src="$rc_sb/Y$sY/M$sM/$base.${sy}_0000z.bin"
  if [[ ! -e "$src" ]]; then
    local py ph; read -r py ph < <("$TICK" "$sy" "$sh" -86400)
    src="$rc_sb/Y${py:0:4}/M${py:4:2}/$base.${py}_0000z.bin"
  fi
  [[ -e "$src" ]] || { echo "(warn) no 00Z catparam near $sy"; return; }

  local ny nh; read -r ny nh < <("$TICK" "$sy" "$sh" 86400)
  local dst="$rc_sb/Y${ny:0:4}/M${ny:4:2}/$base.${ny}_${nh:0:4}z.bin"  # -> *_0000z.bin
  if [[ ! -e "$dst" && ! -L "$dst" ]]; then
    mkdir -p "$(dirname "$dst")"
    ln -s "$src" "$dst"
    echo "Linked catparam: $(basename "$dst") -> $(basename "$src")"
  fi
}

axis_for_tasks() {
  local rdir="${1:-$EXPDIR/run}"
  local nx=""; local ny=""
  if [[ -f "$rdir/LDAS.rc" ]]; then
    nx="$(awk -F: '/^[[:space:]]*NX[[:space:]]*:/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$rdir/LDAS.rc")"
    ny="$(awk -F: '/^[[:space:]]*NY[[:space:]]*:/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$rdir/LDAS.rc")"
  fi
  [[ -z "$nx" ]] && nx=1
  [[ -z "$ny" ]] && ny=1

  if   (( ny > 1 && nx <= 1 )); then echo "NY"
  elif (( nx > 1 && ny <= 1 )); then echo "NX"
  else
    # ambiguous (both 1 or both >1): fall back to EXPDOMAIN
    case "$EXPDOMAIN" in
      CF*) echo "NY" ;;   # CF grids -> NY axis is the 2nd dim
      *)   echo "NX" ;;   # EASE etc.
    esac
    
  fi
}

# Generic NetCDF comparison helper (per file)
#   mode = restart | history | layout
nc_compare() {
  local mode="$1" f1="$2" f2="$3" base="$4"
  local flags

  case "$mode" in
    restart)
      # strict data + metadata + global attrs
      flags="-dmfgMNS"
      ;;
    history)
      # data-only for pass/fail (no metadata)
      flags="-dNM"
      ;;
    layout)
      # tolerant for layout tests (uses NCCMP_FLAGS_TOL)
      flags="$NCCMP_FLAGS_TOL"
      ;;
    *)
      echo "nc_compare: unknown mode '$mode' for $base" >&2
      return 1
      ;;
  esac

  if ! "$NCCMP" $flags "$f1" "$f2"; then
    echo "${mode^^} DIFF: $base"
    return 1
  fi

  return 0
}

#
#
# -------- main ---------------------------------------------------------------
#
echo "== GEOSldas start/stop (model-only) =="

# Read current cap_restart (this is our intended ORIGINAL start)
read START_NYMD START_NHMS < <(cap_restart_start)
echo "Initial start from cap_restart: ${START_NYMD} ${START_NHMS}"

# Keep copies of the original start for the whole script
ORIG_NYMD="$START_NYMD"
# normalize START_NHMS to HHMMSS if needed
if [[ ${#START_NHMS} -eq 4 ]]; then ORIG_NHMS="${START_NHMS}00"; else ORIG_NHMS="$START_NHMS"; fi
REAL_EXPDIR="$EXPDIR"
# Build sandbox and then run entirely inside it
make_sandbox
EXPDIR="$SANDBOX_EXPDIR"         # <---- switch context to sandbox
JOBFILE="$EXPDIR/run/lenkf.j"    # submit the sandbox jobfile

# ---- choose CF (2d/Nx) vs EASE (1d/Nt) by your naming convention ----
# Examples you showed:
#   CF0090x6C_GLOBAL        -> CF
#   SMAP_EASEv2_M36_GLOBAL  -> EASE
if [[ "$EXPDOMAIN" == CF* ]]; then
  GRID_KIND="CF"
else
  GRID_KIND="EASE"
fi

# Install the matching HISTORY template into the sandbox, expects land and landice tile type
if [[ "$GRID_KIND" == "CF" ]]; then
  cp -f "$TEMPLATE_DIR/HISTORY_2d.rc" "$EXPDIR/run/HISTORY.rc"
else
  cp -f "$TEMPLATE_DIR/HISTORY_1d.rc" "$EXPDIR/run/HISTORY.rc"
fi

# If there is no land-ice restart, drop glc collections from HISTORY.rc
if [[ "${HAS_LANDICE:-0}" -eq 0 ]]; then
  echo "Land-only sandbox: removing glc HISTORY collections (no landice_internal_rst)"
  # EASE 1-D glc collection
  sed -i "/tavg06_1d_glc_Nt/d" "$EXPDIR/run/HISTORY.rc" || true
  # CF 2-D glc collection
  sed -i "/tavg06_2d_glc_Nx/d" "$EXPDIR/run/HISTORY.rc" || true
fi

# Ensure LDAS uses this HISTORY.rc
sed -i -E 's|^\s*(HISTORY(_RC)?\s*:\s*).*$|\1HISTORY.rc|' "$EXPDIR/run/LDAS.rc"

# If requested, patch sandbox HISTORY to hourly
if [[ "$HIST_MODE" == "6h" ]]; then
  patch_history_6h
fi

echo "Sanity: checking LDAS.rc for HISTORY settings ..."
grep -nE '^(HISTORY(_RC)?|HISTORY_FILE|DO_HIST|DO_HISTORY)\s*:' "$EXPDIR/run/LDAS.rc" || \
  echo "(note) No explicit HISTORY keys found; GEOSldas may default to run/HISTORY.rc"

echo "HISTORY frequencies (expect 060000):"
grep -nE 'tavg06_(2d_(lfs|lnd|glc)_Nx|1d_(lfs|lnd|glc)_Nt)\.frequency' "$EXPDIR/run/HISTORY.rc" | sed 's/^/  /' || true


# Replace the "finger" command in lenkf.j with a safe user ID
sed -i -E 's|setenv[[:space:]]+MYNAME.*|setenv MYNAME "`id -un`"|' "$EXPDIR/run/lenkf.j"

# Preflight check: verifying sandbox executable and job setup ...
echo "Preflight check: verifying sandbox executable and job setup ..."
ls -l "$EXPDIR/build/bin/GEOSldas.x" "$EXPDIR/build/bin/esma_mpirun" \
  || die "Missing build binaries in sandbox"

# Ensure sandbox starts from the original time and restarts
restore_initial_restarts

# Sanity: is there a catparam at or before the start time?
echo "Checking sandbox rc_out for catparam near $ORIG_NYMD ..."
ls -l "$EXPDIR/output/$EXPDOMAIN/rc_out"/Y*/M*/*"${EXPID}.ldas_catparam."* || true

# T1: single 24h to final time F
# Use GEOS 'tick' for portable datetime math
TICK="$EXPDIR/build/bin/tick"
if [[ ! -x "$TICK" ]]; then
  echo "ERROR: tick not found at $TICK"; exit 1
fi

# +12h (mid) and +24h (final) in seconds
read M_NYMD M_NHMS < <("$TICK" "$START_NYMD" "$START_NHMS" $((12*3600)))
read F_NYMD F_NHMS < <("$TICK" "$START_NYMD" "$START_NHMS" $((24*3600)))

# In patch mode, synthesize 6h CENTER catparam stamps (03/09/15/21Z) in the sandbox
if [[ "$HIST_MODE" == "6h" ]]; then
ensure_catparam_next00z "$START_NYMD" "$START_NHMS"
fi

echo "T1 final: ${F_NYMD} ${F_NHMS}"
set_cap_rc "$T1_HHMMSS" "$F_NYMD" "$F_NHMS"
  echo "T1 final: ${F_NYMD} ${F_NHMS}"
  set_cap_rc "$T1_HHMMSS" "$F_NYMD" "$F_NHMS"

  # ---- Snapshot T1 run dir as template for layout test ----
  T1_TEMPLATE="$EXPDIR/run_T1_template/run"
  if [[ ! -d "$T1_TEMPLATE" ]]; then
    echo "Saving T1 run template to: $T1_TEMPLATE"
    mkdir -p "$(dirname "$T1_TEMPLATE")"
    rsync -a --delete "$EXPDIR/run/" "$T1_TEMPLATE/"
  else
    echo "Using existing T1 run template: $T1_TEMPLATE"
  fi

  submit_and_wait_done
  wait_for_restart_stamp "$F_NYMD" "$F_NHMS"
  collect_set "T1_24h_${F_NYMD}_${F_NHMS}" "$F_NYMD" "$F_NHMS"

if [[ "$HIST_MODE" == "6h" ]]; then
  collect_hist_range "T1_hist_${F_NYMD}_${F_NHMS}" "$START_NYMD" "$START_NHMS" "$F_NYMD" "$F_NHMS"
else
  # daily tavg06: compare at mid (12Z)
  collect_hist_set "T1_hist_${M_NYMD}_${M_NHMS}" "$M_NYMD" "$M_NHMS"
fi
echo "---- HISTORY lines from sandbox log (tail) ----"
# adjust path if your Slurm output names differ
sed -n '1,200p' "$EXPDIR/scratch/GEOSldas_log_txt" | grep -iE 'HISTORY|tavg|CFIO|collection|writing' || true
# HISTORY activity in sandbox log (full scan)
grep -nEi 'HISTORY|tavg06|CFIO|collection|write|wrote' "$EXPDIR/scratch/GEOSldas_log_txt" || echo "  (no HISTORY mentions found)"

# Files anywhere under sandbox
find "$EXPDIR" -maxdepth 6 -type f -name "${EXPID}.tavg06_2d_*_Nx.20*.nc4" -printf '  %p\n' | sort | head -40

# Files in the target cat/ens0000 layout
find "$EXPDIR/output/$EXPDOMAIN/cat/ens0000" -maxdepth 2 -type f -name "${EXPID}.tavg06_2d_*_Nx.20250116_*.nc4" -printf '  %P\n' | sort


# >>> RESET to original start before T2 <<<
restore_initial_restarts

# T2: 12h to mid time M
echo "T2 final (mid): ${M_NYMD} ${M_NHMS}"
set_cap_rc "$T2_HHMMSS" "$M_NYMD" "$M_NHMS"
submit_and_wait_done
wait_for_restart_stamp "$M_NYMD" "$M_NHMS"
collect_set "T2_12h_${M_NYMD}_${M_NHMS}" "$M_NYMD" "$M_NHMS"

if [[ "$HIST_MODE" == "6h" ]]; then
  collect_hist_range "T2_hist_${M_NYMD}_${M_NHMS}" "$START_NYMD" "$START_NHMS" "$M_NYMD" "$M_NHMS"
else
  collect_hist_set "T2_hist_${M_NYMD}_${M_NHMS}" "$M_NYMD" "$M_NHMS"
fi

# Promote T2 restarts as inputs for T3 (continue another 12h)
promote_restarts_from_stamp "$M_NYMD" "$M_NHMS"

# T3: 12h more to same final F
set_cap_rc "$T3_HHMMSS" "$F_NYMD" "$F_NHMS"
submit_and_wait_done
wait_for_restart_stamp "$F_NYMD" "$F_NHMS"
collect_set "T3_12h_from_mid_${F_NYMD}_${F_NHMS}" "$F_NYMD" "$F_NHMS"

if [[ "$HIST_MODE" == "6h" ]]; then
  collect_hist_range "T3_hist_${F_NYMD}_${F_NHMS}" "$M_NYMD" "$M_NHMS" "$F_NYMD" "$F_NHMS"
fi

echo "Compare T1 (24h) vs T3 (12h+12h) at ${F_NYMD}_${F_NHMS}"
if compare_sets "T1_24h_${F_NYMD}_${F_NHMS}" "T3_12h_from_mid_${F_NYMD}_${F_NHMS}"; then
  echo "=========================================."
  echo "start/stop model"
  echo "Test ✅ success"
  echo "RESTARTS: PASS"
  echo "=========================================."
else
  echo "=========================================."
  echo "start/stop model"
  echo "Test ❌ fail"
  echo "RESTARTS: FAIL (differences detected)"
  echo "=========================================."
  exit 2
fi
# --- Sub-daily HISTORY compare (for time-averaged collections) ---
# Supports 6-hour or 1-hour averaging; centers are half the averaging window.
HIST_MODE="${HIST_MODE:-6h}"
case "$HIST_MODE" in
  6h)  HIST_STEP_SEC=21600; HIST_STEP_OFFSET_SEC=10800 ;;  # centers at 03/09/15/21Z (6h means 3h offset)
  1h)  HIST_STEP_SEC=3600;  HIST_STEP_OFFSET_SEC=1800  ;;  # centers at hh:30 (1h means 30min offset)
  *)   echo "(warn) '$HIST_MODE' not supported for restart regression; forcing 6h."
       HIST_MODE="6h"; HIST_STEP_SEC=21600; HIST_STEP_OFFSET_SEC=10800 ;;
esac

# Merge T2 and T3 HISTORY outputs into one stitched 24h set
T23_COMBINED="T23_hist_${F_NYMD}_${F_NHMS}"
mkdir -p "$SETS/$T23_COMBINED"
cp -pn "$SETS/T2_hist_${M_NYMD}_${M_NHMS}"/*.nc4 "$SETS/$T23_COMBINED/" 2>/dev/null || true
cp -pn "$SETS/T3_hist_${F_NYMD}_${F_NHMS}"/*.nc4 "$SETS/$T23_COMBINED/" 2>/dev/null || true

# Sanity check: both sides must contain at least one file
cnt_T1=$(find "$SETS/T1_hist_${F_NYMD}_${F_NHMS}" -maxdepth 1 -type f -name '*.nc4' | wc -l)
cnt_T23=$(find "$SETS/$T23_COMBINED"               -maxdepth 1 -type f -name '*.nc4' | wc -l)
if (( cnt_T1 == 0 || cnt_T23 == 0 )); then
  echo "(error) Empty HISTORY set(s): T1=$cnt_T1 T23=$cnt_T23"; exit 2
fi

echo "Compare HOURLY HISTORY: T1 (24h) vs [T2+T3] (12h+12h)"
if compare_hist_sets "T1_hist_${F_NYMD}_${F_NHMS}" "$T23_COMBINED"; then
  echo "===========================================."
  echo "start/stop model"
  echo "Test ✅ success"
  echo "HISTORY (hourly): PASS"
  echo "===========================================."
else
  echo "===========================================."
  echo "start/stop model"
  echo "Test ❌ fail"
  echo "HISTORY (hourly): FAIL (differences detected)"
  echo "===========================================."
fi

# ---------- T4: layout-invariance (rerun T1 with different ntasks, fully inside sandbox) ----------
# Enable with: RUN_LAYOUT=1 ALT_1D=<int>
#
if [[ "${RUN_LAYOUT:-0}" == "1" ]]; then
  echo "== T4: layout-invariance (T1 rerun with alternate layout, sandbox-only) =="

  : "${ALT_1D:?ALT_1D must be set (e.g., 120)}"

  BASE_EXPDIR="$EXPDIR"           # sandbox root
  BASE_RUN="$BASE_EXPDIR/run"
  BASE_OUT="$BASE_EXPDIR/output/$EXPDOMAIN"
  BASE_RS="$BASE_OUT/rs/ens0000"
  BASE_CAT="$BASE_OUT/cat/ens0000"

  # Prefer a frozen T1 template; required for T4
  T1_TEMPLATE="$BASE_EXPDIR/run_T1_template/run"
  if [[ -d "$T1_TEMPLATE" ]]; then
    CFG_RUN="$T1_TEMPLATE"
    echo "T4: using T1 template run dir: $CFG_RUN"
  else
    echo "T4: ERROR: T1 template not found at $T1_TEMPLATE; skipping layout test."
    return 0
  fi

  # T4 sandbox experiment under sandbox root
  T4_EXPDIR="$BASE_EXPDIR/run_T4"
  T4_RUN="$T4_EXPDIR/run"
  T4_OUT="$T4_EXPDIR/output/$EXPDOMAIN"
  T4_LDAS="$T4_RUN/LDAS.rc"
  mkdir -p "$T4_RUN" "$T4_EXPDIR/scratch" "$T4_OUT" "$T4_OUT/cat/ens0000"
  # Ensure T4 has rc_out copied from sandbox OUT_PATH
  T4_RC_OUT="$T4_OUT/rc_out"
  SB_RC_OUT="$BASE_EXPDIR/output/$EXPDOMAIN/rc_out"

  echo "T4: syncing rc_out from sandbox $SB_RC_OUT -> $T4_RC_OUT"
  mkdir -p "$T4_RC_OUT"
  rsync -a "$SB_RC_OUT/" "$T4_RC_OUT/" || echo "(warn) T4: rc_out rsync failed"


  echo "T4: cloning $CFG_RUN into $T4_RUN"
  rsync -a --delete "$CFG_RUN/" "$T4_RUN/"

  # Copy input/ so T4 has its own input tree; link build/ for binaries only
  mkdir -p "$T4_EXPDIR/input"
  rsync -a --delete "$BASE_EXPDIR/input/" "$T4_EXPDIR/input/"
  ln -sfn "$BASE_EXPDIR/build" "$T4_EXPDIR/build"
 

  # OUT_PATH → T4 OUT (path only)
  sed -i -E "s|^OUT_PATH:[[:space:]]*.*$|OUT_PATH:                         $T4_OUT|" "$T4_RUN/LDAS.rc"
  # Ensure T4 has its own HISTORY.rc identical to T1 template
  cp -f "$CFG_RUN/HISTORY.rc" "$T4_RUN/HISTORY.rc"

  # CAP: point HIST_CF to T4's HISTORY.rc
  if grep -qE '^[[:space:]]*HIST_CF:' "$T4_RUN/CAP.rc"; then
    sed -i -E "s|^[[:space:]]*HIST_CF:[[:space:]]*.*$|HIST_CF: $T4_RUN/HISTORY.rc|" "$T4_RUN/CAP.rc"
  else
    printf 'HIST_CF: %s\n' "$T4_RUN/HISTORY.rc" >> "$T4_RUN/CAP.rc"
  fi

  # T4 jobfile: only change EXPDIR, logs, ntasks/npes (layout only)
  T4_JOB="$T4_RUN/lenkf.j"
  [[ -f "$T4_JOB" ]] || T4_JOB="$(find "$T4_RUN" -maxdepth 1 -type f -name '*.j' | head -n1)"
  [[ -f "$T4_JOB" ]] || die "T4: no jobfile in $T4_RUN"

  # In the jobfile, reset any HISTORY_RC env to point to T4's HISTORY.rc
  sed -i -E '/^[[:space:]]*setenv[[:space:]]+HISTORY_RC/d' "$T4_JOB"
  sed -i -E '/^[[:space:]]*(export[[:space:]]+)?HISTORY_RC=/d' "$T4_JOB"
  echo "setenv HISTORY_RC $T4_RUN/HISTORY.rc" >> "$T4_JOB"

  # EXPDIR → T4_EXPDIR
  sed -i -E "s|^(setenv[[:space:]]+EXPDIR[[:space:]]+).*|\1$T4_EXPDIR|" "$T4_JOB" 2>/dev/null || true
  sed -i -E "s|^(export[[:space:]]+EXPDIR=).*|\1$T4_EXPDIR|"         "$T4_JOB" 2>/dev/null || true

  # Logs → T4 scratch
  sed -i -E "s|^#SBATCH --output=.*|#SBATCH --output=$T4_EXPDIR/scratch/GEOSldas_log_txt|"  "$T4_JOB"
  sed -i -E "s|^#SBATCH --error=.*|#SBATCH --error=$T4_EXPDIR/scratch/GEOSldas_err_txt|"    "$T4_JOB"

  # Slurm ntasks → ALT_1D
  if grep -qE '^#SBATCH[[:space:]]+--ntasks=' "$T4_JOB"; then
    sed -i -E "s|^#SBATCH[[:space:]]+--ntasks=.*|#SBATCH --ntasks=${ALT_1D}|" "$T4_JOB"
  elif grep -qE '^#SBATCH[[:space:]]+-n[[:space:]]+[0-9]+' "$T4_JOB"; then
    sed -i -E "s|^#SBATCH[[:space:]]+-n[[:space:]]+[0-9]+|#SBATCH -n ${ALT_1D}|" "$T4_JOB"
  else
    sed -i '1a #SBATCH --ntasks='"$ALT_1D" "$T4_JOB"
  fi

  # Optional resource tweaks (can be removed if you prefer)
  grep -qE '^#SBATCH[[:space:]]+--constraint=' "$T4_JOB" || sed -i '1a #SBATCH --constraint=mil' "$T4_JOB"
  if (( ALT_1D <= 126 )); then
    grep -qE '^#SBATCH[[:space:]]+--nodes=' "$T4_JOB" || sed -i '1a #SBATCH --nodes=1' "$T4_JOB"
  fi

  # MPI launch mirrors for npes
  sed -i -E "s/(--npes_model[[:space:]]+)[0-9]+/\1${ALT_1D}/g" "$T4_JOB"
  sed -i -E "/(mpirun|esma_mpirun)/ s/(-np[[:space:]]+)[0-9]+/\1${ALT_1D}/" "$T4_JOB"
  sed -i -E 's/\$RUN_CMD[[:space:]]+\$?total_npes([[:space:]]+)/\$RUN_CMD '"$ALT_1D"'\1/' "$T4_JOB"
  sed -i -E 's/^[[:space:]]*set[[:space:]]+total_npes[[:space:]]*=.*/set total_npes = '"$ALT_1D"'/' "$T4_JOB"
  grep -qE '(^|[[:space:]])setenv[[:space:]]+GEOS_NP' "$T4_JOB" \
    && sed -i -E "s/^[[:space:]]*setenv[[:space:]]+GEOS_NP[[:space:]]+.*/setenv GEOS_NP ${ALT_1D}/" "$T4_JOB" \
    || printf 'setenv GEOS_NP %s\n' "$ALT_1D" >> "$T4_JOB"
  grep -qE '(^|[[:space:]])setenv[[:space:]]+NTASKS_?MODEL' "$T4_JOB" \
    && sed -i -E "s/^[[:space:]]*setenv[[:space:]]+NTASKS_?MODEL[[:space:]]+.*/setenv NTASKS_MODEL ${ALT_1D}/" "$T4_JOB" \
    || printf 'setenv NTASKS_MODEL %s\n' "$ALT_1D" >> "$T4_JOB"

  # Remove any IM/JM overrides (dims unchanged)
  sed -i -E '/(^|[[:space:]])(setenv[[:space:]]+|export[[:space:]]+)?IM_WORLD([[:space:]]|=)/d' "$T4_JOB"
  sed -i -E '/(^|[[:space:]]+)(setenv[[:space:]]+|export[[:space:]]+)?JM_WORLD([[:space:]]|=)/d' "$T4_JOB"

  # 3) Axis & layout file (CF → NY/JMS.rc, EASE → NX/IMS.rc)
  AXIS_KEY="$(axis_for_tasks "$BASE_RUN")"   # "NX" or "NY"
  if [[ "$AXIS_KEY" == "NY" ]]; then
    T4_LAYOUT_FILE="$T4_RUN/JMS.rc"
  else
    T4_LAYOUT_FILE="$T4_RUN/IMS.rc"
  fi

  # --- Prebuild layout (JMS/IMS) so the model can start with ALT_1D ---
  # TILE_TYPES: read from LDAS.rc only; don't change it for now
  TILE_TYPES="$(
    awk -F: '/^[[:space:]]*TILE_TYPES[[:space:]]*:/{
      v=$2; gsub(/^[ \t]+|[ \t]+$/,"",v); gsub(/[ \t]+/,"_",v); print v; exit
    }' "$T4_LDAS"
  )"
  
  if [[ -z "$TILE_TYPES" || "$TILE_TYPES" == "_" ]]; then
    echo "❌ TILE_TYPES not found in $T4_LDAS; add TILE_TYPES: to LDAS.rc or re-enable auto-detect."
    exit 2
  fi
  
  LOG="$T4_RUN/preprocess_ldas.log"
  PPL="$T4_EXPDIR/build/bin/preprocess_ldas.x"
  G5C="$BASE_EXPDIR/build/bin/g5_modules"   # or REAL_EXPDIR, whichever you use elsewhere
  OPT_FILE="$T4_RUN/LDAS.optimized"
  
  if [[ -x "$PPL" && -f "$G5C" ]]; then
    echo "PREPROC: pre-building layout for ALT=$ALT_1D (TILE_TYPES=$TILE_TYPES) → $LOG"
    csh -f -c "source $G5C; cd $T4_RUN; \
      $PPL optimize \
        $T4_EXPDIR/input/tile.data \
        $ALT_1D \
        $OPT_FILE \
        $T4_RUN \
        $TILE_TYPES" |& tee "$LOG"
  else
    echo "❌ Missing preprocess or env: PPL=$PPL  g5_modules=$G5C"
    exit 2
  fi
  
  # Verify that the generated layout file header matches ALT_1D
  hdr="$(awk 'NR==1{print $1}' "$T4_LAYOUT_FILE" 2>/dev/null || echo '')"
  [[ "$hdr" == "$ALT_1D" ]] || { echo "❌ $(basename "$T4_LAYOUT_FILE") header=$hdr, expected $ALT_1D"; exit 2; }
  # --- end prebuild ---

  # LDAS.NX/NY line → ALT_1D (active axis only)
  sed -i -E "s~^([[:space:]]*$AXIS_KEY[[:space:]]*:[[:space:]]*).*~\1 ${ALT_1D}~" "$T4_LDAS"

  # Restarts: copy FROM sandbox rs (the ones T1 started from) into private T4 rs tree
  ORIG_Y="${ORIG_NYMD:0:4}"; ORIG_M="${ORIG_NYMD:4:2}"
  T4_RS_ROOT="$T4_OUT/rs/ens0000"
  T4_RS_DIR="$T4_RS_ROOT/Y$ORIG_Y/M$ORIG_M"
  mkdir -p "$T4_RS_DIR" "$T4_EXPDIR/input/restart"

  echo "T4: copying sandbox start restarts from $BASE_RS to $T4_RS_DIR for ${ORIG_NYMD}_${ORIG_NHMS}"
  
  for comp in catch landice; do
    # Land-only case: skip landice restart entirely
    if [[ "$comp" == "landice" && "${HAS_LANDICE:-0}" -eq 0 ]]; then
      echo "  T4: skipping landice restart (HAS_LANDICE=0)"
      continue
    fi
  
    base="$BASE_RS/Y$ORIG_Y/M$ORIG_M/$EXPID.${comp}_internal_rst"
    src=""
  
    [[ -f "${base}.${ORIG_NYMD}_${ORIG_NHMS}"           ]] && src="${base}.${ORIG_NYMD}_${ORIG_NHMS}"
    [[ -z "$src" && -f "${base}.${ORIG_NYMD}_${ORIG_NHMS:0:4}" ]] && src="${base}.${ORIG_NYMD}_${ORIG_NHMS:0:4}"
  
    if [[ -z "$src" ]]; then
      echo "❌ T4: missing sandbox $comp restart ${ORIG_NYMD}_${ORIG_NHMS} at $base.*"
      # For catch, this is fatal; landice should never hit here if HAS_LANDICE=0
      exit 2
    fi
  
    tgt="$T4_RS_DIR/$(basename "$src")"
    cp -p "$src" "$tgt"
    ln -rsf "$tgt" "$T4_EXPDIR/input/restart/${comp}_internal_rst"
    echo "  -> copied $comp restart to T4: $tgt"
  done

  VEG_SRC="$BASE_RS/$EXPID.vegdyn_internal_rst"
  if [[ -f "$VEG_SRC" ]]; then
    VEG_TGT="$T4_RS_ROOT/$(basename "$VEG_SRC")"
    cp -p "$VEG_SRC" "$VEG_TGT"
    ln -rsf "$VEG_TGT" "$T4_EXPDIR/input/restart/vegdyn_internal_rst"
    echo "  -> copied vegdyn_internal_rst to T4: $VEG_TGT"
  else
    echo "(warn) T4: vegdyn_internal_rst not found at $VEG_SRC"
  fi
  # CAP.rc & cap_restart: copy exactly from T1 template (test window & segments)
  cp -p "$CFG_RUN/CAP.rc"      "$T4_RUN/CAP.rc"
  cp -p "$CFG_RUN/cap_restart" "$T4_RUN/cap_restart" 2>/dev/null || true

  echo "T4 CAP.rc window:"
  grep -E '^(BEG_DATE|END_DATE|JOB_SGMT):' "$T4_RUN/CAP.rc" || true

  # Run T4 inside sandbox
  echo ">>> T4 launch (layout patched to ${ALT_1D}):"
  grep -nE 'mpirun|esma_mpirun|--npes_model| -np |GEOS_NP|NTASKS_MODEL|total_npes' "$T4_JOB" | sed 's/^/  /' || true

  SAVE_EXPDIR="$EXPDIR"
  EXPDIR="$T4_EXPDIR"
  JOBFILE="$T4_JOB"
  submit_and_wait_done
  EXPDIR="$SAVE_EXPDIR"

  # HISTORY & restart compare: T1 (sandbox) vs T4 (sandbox)
  TICK="$BASE_EXPDIR/build/bin/tick"
  if [[ ! -x "$TICK" ]]; then
    echo "(warn) tick not found at $TICK; skipping T4 compares"
  else
    read BASE_BEG_Y BASE_BEG_H < <("$TICK" "$ORIG_NYMD" "$ORIG_NHMS" -86400)

    _has_hist() {
      local cat="$1" y="$2" patt="$3"
      find "$cat" -type f -name "${EXPID}.*.${y}_${patt}z*.nc4" -print -quit 2>/dev/null | grep -q .
    }

    CADENCE="hourly"
    _has_hist "$BASE_CAT" "$BASE_BEG_Y" "1200" && CADENCE="daily12Z"
    if _has_hist "$BASE_CAT" "$BASE_BEG_Y" "03??" || _has_hist "$BASE_CAT" "$BASE_BEG_Y" "0900" \
       || _has_hist "$BASE_CAT" "$BASE_BEG_Y" "1500" || _has_hist "$BASE_CAT" "$BASE_BEG_Y" "2100"; then
      CADENCE="6h"
    fi
    echo "Baseline (T1) HISTORY cadence: $CADENCE"

    BASE_TAG="T1_layout_hist_${F_NYMD}_${F_NHMS}"
    T4_TAG="T4_layout_hist_${F_NYMD}_${F_NHMS}"

    # collect baseline (T1) and T4 HISTORY from sandbox
    if [[ "$CADENCE" == "daily12Z" ]]; then
      read YY12 _ < <("$TICK" "$F_NYMD" "$F_NHMS" -43200)
      HH12="120000"
      SAVE="$EXPDIR"; EXPDIR="$BASE_EXPDIR"; collect_hist_set "$BASE_TAG" "$YY12" "$HH12"; EXPDIR="$SAVE"
      SAVE="$EXPDIR"; EXPDIR="$T4_EXPDIR";   collect_hist_set "$T4_TAG"   "$YY12" "$HH12"; EXPDIR="$SAVE"
    else
      SAVE="$EXPDIR"; EXPDIR="$BASE_EXPDIR"; collect_hist_range "$BASE_TAG" "$BASE_BEG_Y" "$BASE_BEG_H" "$F_NYMD" "$F_NHMS"; EXPDIR="$SAVE"
      SAVE="$EXPDIR"; EXPDIR="$T4_EXPDIR";   collect_hist_range "$T4_TAG"   "$BASE_BEG_Y" "$BASE_BEG_H" "$F_NYMD" "$F_NHMS"; EXPDIR="$SAVE"
    fi

    echo "→ HISTORY (tolerant) T1 vs T4"
    if ! compare_hist_sets "$BASE_TAG" "$T4_TAG"; then
      echo "❌ HISTORY differs under alternate layout"; exit 2
    fi
    echo "✔ HISTORY tolerant compare OK"

    # Restart compare at final time (strict)
    SAVE="$EXPDIR"; EXPDIR="$BASE_EXPDIR"; collect_set "T1_layout_24h_${F_NYMD}_${F_NHMS}" "$F_NYMD" "$F_NHMS"; EXPDIR="$SAVE"
    SAVE="$EXPDIR"; EXPDIR="$T4_EXPDIR";   collect_set "T4_layout_24h_${F_NYMD}_${F_NHMS}" "$F_NYMD" "$F_NHMS";   EXPDIR="$SAVE"

    echo "→ Restarts (strict) T1 vs T4"
    if ! compare_sets "T1_layout_24h_${F_NYMD}_${F_NHMS}" "T4_layout_24h_${F_NYMD}_${F_NHMS}"; then
      echo "❌ Restarts differ under layout"; exit 2
    fi
    echo "✔ Restarts identical (layout-invariance)"

    # Final layout report
    get_axis_val() {
      local file="$1" axis="$2"
      awk -F: -v k="$axis" '
        $1 ~ ("^[[:space:]]*" k "[[:space:]]*$") { gsub(/[[:space:]]/,"",$2); print $2; exit }' "$file"
    }
    b_nx="$(get_axis_val "$BASE_RUN/LDAS.rc" NX)"
    b_ny="$(get_axis_val "$BASE_RUN/LDAS.rc" NY)"
    t_nx="$(get_axis_val "$T4_LDAS" NX)"
    t_ny="$(get_axis_val "$T4_LDAS" NY)"
    layout_file="$(basename "$T4_LAYOUT_FILE" 2>/dev/null || echo unknown)"

    echo "=============================================="
    echo "T4 layout test ✅ success"
    echo "T1 grid:   NX=$b_nx  NY=$b_ny"
    echo "T4 grid:   NX=$t_nx  NY=$t_ny"
    echo "Axis:      $AXIS_KEY  (changed → $ALT_1D)"
    echo "Layout:    $layout_file (first integer = $ALT_1D)"
    echo "Changed run files: LDAS.rc, lenkf.j, $layout_file"
    echo "Restarts (24h) & HISTORY (24h) match T1"
    echo "=============================================="
  fi
fi
# ---------- end T4: layout-invariance ----------

echo "Done."

# Comment to see sandbox after a pass
 rm -rf "$SANDBOX_ROOT"     # or: rm -rf "$REAL_EXPDIR/regress/sandbox"
 echo "Cleaning up sandbox: $SANDBOX_ROOT"

