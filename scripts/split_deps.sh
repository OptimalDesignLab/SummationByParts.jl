#!/usr/bin/env bash
set -euo pipefail

# usage: ./scripts/split_deps.sh --dry-run 

# Ensure running in Bash
if [ -z "$BASH_VERSION" ]; then
    exec bash "$0" "$@"
fi

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then DRY_RUN=1; fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
# if [[ "$CURRENT_BRANCH" == "main" || "$CURRENT_BRANCH" == "master" ]]; then
#     echo "Refusing to run on branch '$CURRENT_BRANCH'. Switch to a feature branch." >&2
#     exit 1
# fi

SCAN_DIRS=(src examples test docs)

# Extract deps from Project.toml
DEPS=()
while IFS= read -r line; do
    DEPS+=("$line")
done < <(awk '/^\[deps\]/{flag=1;next}/^\[/{flag=0}flag && NF{gsub(/ *=.*$/, "", $0); print $0}' Project.toml | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

[[ ${#DEPS[@]} -eq 0 ]] && { echo "No dependencies found in Project.toml [deps]. Exiting."; exit 0; }

# Initialize usage tracking
USAGES=()
for pkg in "${DEPS[@]}"; do
    USAGES+=("$pkg=")
done

normalize_token() {
    sed -E 's/[:,;].*$//; s/^[[:space:]]*//; s/[[:space:]]*$//; s/^\.//'
}

scan_file() {
    local file="$1"
    while IFS= read -r line; do
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        [[ "$line" =~ using|import ]] || continue
        line=$(echo "$line" | sed -E 's/#.*$//; s/^[[:space:]]*(using|import)[[:space:]]+//')
        IFS=',:' read -ra tokens <<< "$line"
        for token in "${tokens[@]}"; do
            token=$(echo "$token" | normalize_token)
            key_token=$(printf '%s' "$token" | tr '[:upper:]' '[:lower:]' | sed 's/_//g')
            for pkg in "${DEPS[@]}"; do
                key_pkg=$(printf '%s' "$pkg" | tr '[:upper:]' '[:lower:]' | sed 's/_//g')
                if [[ "$key_token" == "$key_pkg" ]]; then
                    for i in "${!USAGES[@]}"; do
                        key=${USAGES[$i]%%=*}
                        value=${USAGES[$i]#*=}
                        if [[ "$key" == "$pkg" ]]; then
                            [[ -z "$value" ]] && USAGES[$i]="$pkg=$file" || [[ ! ",$value," =~ ",$file," ]] && USAGES[$i]="$pkg=$value,$file"
                            break
                        fi
                    done
                fi
            done
        done
    done < "$file"
}

scan_notebook() {
    local file="$1"
    if ! command -v jq >/dev/null 2>&1; then
        echo "Warning: jq not installed, skipping notebook $file"
        return
    fi
    jq -e . "$file" >/dev/null 2>&1 || { echo "Warning: $file is not valid JSON, skipping"; return; }

    local tmpfile
    tmpfile=$(mktemp)
    jq -r '.cells[] | select(.cell_type=="code") | .source[]' "$file" > "$tmpfile"

    while IFS= read -r line; do
        line=$(echo "$line" | sed -E 's/^"//;s/"$//;s/\\n$//')
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        [[ "$line" =~ using|import ]] || continue
        line=$(echo "$line" | sed -E 's/#.*$//')
        for token in $line; do
            token=$(echo "$token" | normalize_token)
            key_token=$(printf '%s' "$token" | tr '[:upper:]' '[:lower:]' | sed 's/_//g')
            for i in "${!DEPS[@]}"; do
                key_pkg=$(printf '%s' "${DEPS[$i]}" | tr '[:upper:]' '[:lower:]' | sed 's/_//g')
                if [[ "$key_token" == "$key_pkg" ]]; then
                    key=${DEPS[$i]}
                    value=${USAGES[$i]#*=}
                    [[ -z "$value" ]] && USAGES[$i]="$key=$file" || [[ ! ",$value," =~ ",$file," ]] && USAGES[$i]="$key=$value,$file"
                fi
            done
        done
    done < "$tmpfile"
    rm -f "$tmpfile"
}

# Scan all directories
for d in "${SCAN_DIRS[@]}"; do
    [[ ! -d "$d" ]] && continue
    while IFS= read -r file; do
        ext="${file##*.}"
        case "$ext" in
            jl) scan_file "$file" ;;
            ipynb) scan_notebook "$file" ;;
        esac
    done < <(find "$d" -type f \( -name '*.jl' -o -name '*.ipynb' \) 2>/dev/null)
done

# Determine packages
MOVE_PKGS=()
PKG_TARGETS=()
CORE_PKGS=()
UNMOVED_PKGS=()
for i in "${!USAGES[@]}"; do
    pkg=${USAGES[$i]%%=*}
    usage=${USAGES[$i]#*=}
    if [[ "$usage" == *"src"* ]]; then
        CORE_PKGS+=("$pkg")
    elif [[ -n "$usage" ]]; then
        target=""
        for t in gpu_api examples io eqns test docs; do
            [[ "$usage" == *"$t"* ]] && target="$t" && break
        done
        [[ -z "$target" ]] && target="examples"
        MOVE_PKGS+=("$pkg")
        PKG_TARGETS+=("$pkg=$target")
    else
        UNMOVED_PKGS+=("$pkg")
    fi
done

# Display summary
echo "Detected ${#DEPS[@]} deps in Project.toml."
echo "Core packages (used in src) [${#CORE_PKGS[@]}]:"
for pkg in "${CORE_PKGS[@]}"; do echo "  - $pkg"; done
echo "Packages to move out of root [${#MOVE_PKGS[@]}]:"
for i in "${!MOVE_PKGS[@]}"; do
    pkg=${MOVE_PKGS[$i]}
    target=${PKG_TARGETS[$i]#*=}
    echo "  - $pkg -> target='$target'"
done
echo "Deps in Project.toml but not in src and not moved [${#UNMOVED_PKGS[@]}]:"
if [[ ${#UNMOVED_PKGS[@]} -gt 0 ]]; then
    for pkg in "${UNMOVED_PKGS[@]}"; do
        echo "  - $pkg"
    done
fi

[[ ${#MOVE_PKGS[@]} -eq 0 ]] && { echo "Nothing to move. Exiting."; exit 0; }
[[ $DRY_RUN -eq 1 ]] && { echo "Dry run complete. No changes made."; exit 0; }

# Build Julia script
JULIA_SCRIPT=$(mktemp)
cat > "$JULIA_SCRIPT" <<'JULIA'
using Pkg
Pkg.activate(".")
JULIA

echo -n "Pkg.rm([" >> "$JULIA_SCRIPT"
for i in "${!MOVE_PKGS[@]}"; do
    pkg="${MOVE_PKGS[$i]}"
    if [[ $i -gt 0 ]]; then
        echo -n "," >> "$JULIA_SCRIPT"
    fi
    printf '"%s"' "$pkg" >> "$JULIA_SCRIPT"
done
echo "]);" >> "$JULIA_SCRIPT"

for i in "${!MOVE_PKGS[@]}"; do
    pkg=${MOVE_PKGS[$i]}
    target=${PKG_TARGETS[$i]#*=}
    echo "Pkg.activate(\"$target\")" >> "$JULIA_SCRIPT"
    echo "Pkg.add(\"$pkg\")" >> "$JULIA_SCRIPT"
done
echo "Pkg.instantiate();" >> "$JULIA_SCRIPT"

# Execute Julia script
julia "$JULIA_SCRIPT"
rm -f "$JULIA_SCRIPT"

# Commit
# git add Project.toml Manifest.toml
# git commit -m "chore(env): move detected non-core deps to targets"
# git push -u origin "$CURRENT_BRANCH"

echo "Done: moved packages and pushed branch '$CURRENT_BRANCH'."
