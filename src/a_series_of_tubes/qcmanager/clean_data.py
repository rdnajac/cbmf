import csv


def read_csv(filename):
    with open(filename, newline="") as file:
        reader = csv.DictReader(file)
        return list(reader), reader.fieldnames


def clean_data(data, headers):
    columns_to_remove = {"File type", "Encoding"}
    columns_to_move_to_end = ["Total Deduplicated Percentage"]

    unique_values = {col: set() for col in columns_to_remove}
    for row in data:
        for col in columns_to_remove:
            unique_values[col].add(row[col])
            if len(unique_values[col]) > 1:
                unique_values[col] = None

    columns_to_remove = {
        col for col, values in unique_values.items() if values is not None
    }
    remaining_headers = [col for col in headers if col not in columns_to_remove]

    for col in columns_to_move_to_end:
        if col in remaining_headers:
            remaining_headers.remove(col)
            remaining_headers.append(col)

    cleaned_data = [{k: row[k] for k in remaining_headers} for row in data]

    return cleaned_data, remaining_headers


def merge_data(clean_data, clean_headers, flagstat_data, flagstat_headers):
    flagstat_map = {row["Sample ID"]: row for row in flagstat_data}
    extra_headers = [h for h in flagstat_headers if h not in {"Sample ID"}]

    for row in clean_data:
        flagstat_row = flagstat_map.get(row["Sample ID"])
        if flagstat_row:
            for header in extra_headers:
                row[header] = flagstat_row.get(header, "")

    return clean_data, clean_headers + extra_headers


def write_csv(data, headers, filename):
    with open(filename, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=headers)
        writer.writeheader()
        writer.writerows(data)


def main():
    input_csv = "performance_metrics.csv"
    flagstat_csv = "flagstat.csv"
    output_csv = "final_report.csv"

    data, headers = read_csv(input_csv)
    cleaned_data, cleaned_headers = clean_data(data, headers)

    flagstat_data, flagstat_headers = read_csv(flagstat_csv)
    # print(f"len(flagstat_data): {len(flagstat_data)}")
    # print(f"len(flagstat_headers): {len(flagstat_headers)}")
    # print(f"flagstat_headers: {flagstat_headers}")
    final_data, final_headers = merge_data(
        cleaned_data, cleaned_headers, flagstat_data, flagstat_headers
    )
    write_csv(final_data, final_headers, output_csv)


if __name__ == "__main__":
    main()
