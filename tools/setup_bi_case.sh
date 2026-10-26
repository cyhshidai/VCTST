#!/bin/bash
# Build a ready-to-run ABACUS vc-NEB case from the two common-cell POSCARs
# produced by tools/build_common_cell.py.
#
#   tools/setup_bi_case.sh <case_dir> <build_dir> [np]
#
#   case_dir   where the case is created (default ./bi_case)
#   build_dir  output of tools/build_common_cell.py (default ./build)
#   np         MPI ranks for `command:` (default 14)
#
# It only writes files; it does not run anything.  Next steps:
#   cd <case_dir> && vcneb doctor neb.yaml && vcneb run neb.yaml --dry-run
set -euo pipefail

CASE=${1:-./bi_case}
BUILD=${2:-./build}
NP=${3:-14}

mkdir -p "$CASE"
cp "$BUILD/Bi_A7_common.vasp" "$CASE/Bi_A7c.vasp"
cp "$BUILD/Bi_II_common.vasp" "$CASE/Bi_IIc.vasp"

cat > "$CASE/neb.yaml" <<EOF
# Bi A7 (R-3m) -> Bi-II (C2/m) at 2.55 GPa, both ends in one 12-atom cell.
project:
  name: Bi_A7_to_BiII_2p55GPa
  workdir: ./run

structure:
  initial: Bi_A7c.vasp
  final: Bi_IIc.vasp

calc:
  engine: abacus
  command: mpirun -np ${NP} abacus
  mpi: ${NP}
  pseudo_dir: .                  # UPF pseudopotentials next to this file
  data_dir: .                    # .orb orbitals
  timeout: 21600
  pseudopotentials:
    Bi: Bi_ONCV_PBE-1.0.upf
  orbitals:
    Bi: Bi_gga_10au_100Ry_2s2p2d1f.orb
  kpts: {density: 0.03}          # ONE frozen mesh for the whole path
  params:
    ecutwfc: 100
    basis_type: lcao
    ks_solver: genelpa
    dft_functional: pbe
    smearing_method: mp
    smearing_sigma: 0.01
    scf_thr: 1.0e-7
    scf_nmax: 100
    symmetry: 0

neb:
  images: 7
  interpolation: idpp
  align_cell: true
  min_distance: 2.0
  prune: auto
  variable_cell: true
  pressure: 2.55                 # GPa -> H = E + P*V; keep PRESS at 0 in ABACUS
  spring_k: 0.1
  f_th: 0.05
  max_steps: 400
  dt: 0.1
  max_step: 0.2
  cell_max_step: 0.02
  climbing_image: true
  relax_endpoints: true
  endpoint_fmax: 0.05
  endpoint_steps: 300
  parallel_images: 1

output:
  plot: true
EOF

echo "wrote $CASE/neb.yaml"
echo "now copy the UPF/orbital files into $CASE, then:"
echo "  cd $CASE && vcneb doctor neb.yaml && vcneb run neb.yaml --dry-run"
