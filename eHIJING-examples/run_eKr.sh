#!/bin/bash
set -euo pipefail

# Define simulation arguments
Nevents=100
K=4.0
M=1 # Generalized HT:1,  Higher-Twist:0, both in the soft gluon emission limit.
Z=36
A=84

# Define paths
executable_path="/SMASH/mnt/src/ehijing/eHIJING-examples/build/ehijing-space-time-position"
config_path="/SMASH/mnt/src/ehijing/eHIJING-examples/HERMES.setting"
output_dir="/SMASH/mnt/run/ehijing/events"
table_path="/SMASH/mnt/src/ehijing/eHIJING-examples/Tables/${K}"

# Create output directories
mkdir -p "${output_dir}"
mkdir -p "${table_path}"

# Quick sanity checks
[[ -x "${executable_path}" ]] || { echo "ERROR: not executable: ${executable_path}" >&2; exit 1; }
[[ -f "${config_path}" ]]      || { echo "ERROR: missing config: ${config_path}" >&2; exit 1; }

# Run executable
"${executable_path}" "${Nevents}" "${Z}" "${A}" "${M}" "${K}" "${table_path}" "${output_dir}" "${config_path}" # > /dev/null 2>&1



