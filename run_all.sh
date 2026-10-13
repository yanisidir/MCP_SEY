#!/bin/bash
#
# Lance l'ensemble des balayages de l'etude, puis les analyses correspondantes.
#
#   ./run_all.sh            simulations puis analyses
#   ./run_all.sh --sim      simulations seules
#   ./run_all.sh --analyse  analyses seules, sur les .root deja presents
#
# Compter plusieurs heures : le cout d'un point croit avec le gain.

set -euo pipefail

# Se placer dans le dossier du script : les macros utilisent des chemins
# relatifs a la racine du projet.
cd "$(dirname "$0")"

EXEC="./build/MCP_SEY"
MACROS="./macros"
LOGS="./logs"

mode="${1:---tout}"

if [[ "$mode" != "--tout" && "$mode" != "--sim" && "$mode" != "--analyse" ]]; then
    echo "Usage : $0 [--tout | --sim | --analyse]" >&2
    exit 2
fi

mkdir -p "$LOGS" root_files figures

# Une etape qui echoue n'interrompt pas les autres : elles sont independantes.
# Les echecs sont recapitules a la fin, et le code de retour les reflete.
echecs=()

# ---------------------------------------------------------------- simulation

Simuler() {
    local titre="$1" macro="$2" journal="$3"

    if [[ ! -f "$MACROS/$macro" ]]; then
        echo "  macro introuvable : $MACROS/$macro" >&2
        return 1
    fi

    printf '%-42s ' "===== $titre"
    local debut=$SECONDS
    if ! "$EXEC" "$MACROS/$macro" > "$LOGS/$journal" 2>&1; then
        echo "ECHEC"
        echo "--- dernieres lignes de $LOGS/$journal ---" >&2
        tail -15 "$LOGS/$journal" >&2
        return 1
    fi
    printf 'termine en %d min %02d s\n' $(( (SECONDS-debut)/60 )) $(( (SECONDS-debut)%60 ))
}

if [[ "$mode" != "--analyse" ]]; then
    if [[ ! -x "$EXEC" ]]; then
        echo "Executable absent. Compiler d'abord :" >&2
        echo "  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j\$(nproc)" >&2
        exit 1
    fi

    Simuler "Tension"                    scan_voltage.mac     voltage.log     || echecs+=("simulation tension")
    Simuler "Epaisseur"                  scan_thickness.mac   thickness.log   || echecs+=("simulation epaisseur")
    Simuler "Epaisseur a champ constant" scan_field.mac       field.log       || echecs+=("simulation champ constant")
    Simuler "Geometrie de Wu"            scan_wu.mac          wu.log          || echecs+=("simulation Wu")
fi

# ---------------------------------------------------------------- analyse

Analyser() {
    local titre="$1" appel="$2" journal="$3"

    printf '%-42s ' "----- $titre"
    if ! root -l -b -q "$appel" > "$LOGS/$journal" 2>&1; then
        echo "ECHEC"
        tail -10 "$LOGS/$journal" >&2
        return 1
    fi
    echo "figures ecrites"
}

if [[ "$mode" != "--sim" ]]; then
    Analyser "Tension" \
        'analysis/scan_voltage.C("voltage_*200events.root","figures/scan_voltage")' \
        a_voltage.log || echecs+=("analyse tension")
    Analyser "Epaisseur" \
        'analysis/scan_thickness.C("root_files/scanThickness_*.root","figures/scan_thickness")' \
        a_thickness.log || echecs+=("analyse epaisseur")
    Analyser "Epaisseur a champ constant" \
        'analysis/scan_thickness.C("root_files/scanField_*.root","figures/scan_field")' \
        a_field.log || echecs+=("analyse champ constant")
    Analyser "Geometrie de Wu" \
        'analysis/gain_voltage.C("root_files/scanWu_*.root","figures/gain_voltage_wu")' \
        a_wu.log || echecs+=("analyse Wu")
fi

if [[ ${#echecs[@]} -eq 0 ]]; then
    echo "===== Termine sans erreur. Journaux dans $LOGS/, figures dans figures/ ====="
else
    echo "===== Termine avec ${#echecs[@]} echec(s) ====="
    printf '  %s\n' "${echecs[@]}"
    echo "Voir les journaux dans $LOGS/."
    exit 1
fi
