micromamba run -n "$ENV_NAME" fastqc -t "$(nproc)" --memory "$MEMSIZE" --noextract -o "$out_dir" "$input_dir"/*

