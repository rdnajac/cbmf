#!/usr/bin/env sh

err() { printf "\033[1;31mError: %s\033[0m\n" "$1" >&2; }

err_exit() {
	err "$1"
	exit "${2:-1}"
}

SUCCESS() {
	if [ "$(tput colors)" -eq 256 ]; then
		printf '🌈\033[38;5;196mS\033[38;5;208mU\033[38;5;226mC\033[38;5;46mC\033[38;5;21mE\033[38;5;93mS\033[38;5;163mS\033[0m!✨\n'
	elif [ "$(tput colors)" -ge 8 ]; then
		printf '🌈\033[31mS\033[33mU\033[32mC\033[36mC\033[34mE\033[35mS\033[31mS\033[0m!✨\n'
	else
		printf 'SUCCESS!\n'
	fi
}

detect_platform_architecture() {
	case "$(uname)-$(uname -m)" in
	Linux-x86_64) echo "linux-64" ;;
	Linux-aarch64) echo "linux-aarch64" ;;
	Linux-ppc64le) echo "linux-ppc64le" ;;
	Darwin-x86_64) echo "osx-64" ;;
	Darwin-arm64) echo "osx-arm64" ;;
	*) err_exit "Unsupported platform-architecture combination" ;;
	esac
}

detect_shell() {
	DETECTED_SHELL=$(ps -p $$ -o comm=)
	DETECTED_SHELL=${DETECTED_SHELL#-} # Remove leading dash if present
	case "$DETECTED_SHELL" in
	bash | zsh | fish | xonsh) echo "$DETECTED_SHELL" ;;
	*) echo "$SHELL" | xargs basename ;;
	esac
}

initialize_micromamba_shell() {
	MAMBA_ROOT_PREFIX="$1"
	MAMBA_EXE="$2"
	DETECTED_SHELL=$(detect_shell)

	"$MAMBA_EXE" shell init -s "$DETECTED_SHELL" -p "$MAMBA_ROOT_PREFIX" || err_exit "Failed to initialize micromamba shell"

	case "$DETECTED_SHELL" in
	bash) RCFILE="$HOME/.bashrc" ;;
	zsh) RCFILE="$HOME/.zshrc" ;;
	fish) RCFILE="$HOME/.config/fish/config.fish" ;;
	*) err_exit "Unsupported shell for automatic initialization: $DETECTED_SHELL" ;;
	esac

	echo "Micromamba shell initialized. Please run 'source $RCFILE' or restart your terminal to apply changes."
}

write_condarc() {
	cat << EOF > ~/.condarc
channels:
  - bioconda
  - conda-forge
  - nodefaults
channel_priority: strict
EOF
}

cbmf_default_env() {
$MAMBA_EXE create -n cbmf -c bioconda -c conda-forge \
fastqc=0.11.9 hisat2=2.2.1 bwa=0.7.17 bowtie2=2.4.2 \
samtools=1.13 htslib=1.13 bcftools=1.13 \
stringtie=2.1.7 subread=2.0.1
}

