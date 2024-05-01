#! /usr/bin/env python3
import csv
import sys

# Columns to move to the end
move_to_end = ["Total Deduplicated Percentage"]

# Columns to remove if they contain the same value throughout
remove_if_same = ["File type", "Encoding"]


def clean_data(csv_file, output_file):
    """Clean the data in the CSV file."""
    with open(csv_file, "r", newline="") as f:
        reader = csv.reader(f)
        headers = next(reader)

        columns_to_remove = set()
        columns_data = {header: [] for header in headers}

        for row in reader:
            for h, value in zip(headers, row):
                columns_data[h].append(value)

        for col in remove_if_same:
            if all(value == columns_data[col][0] for value in columns_data[col]):
                print(
                    f'Removing  {col} because all entries are "{columns_data[col][0]}"'
                )
                columns_to_remove.add(col)

        new_headers = [col for col in headers if col not in columns_to_remove]
        for col in move_to_end:
            if col in new_headers:
                new_headers.remove(col)
                new_headers.append(col)

        with open(output_file, "w", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(new_headers)

            for row in zip(*[columns_data[col] for col in new_headers]):
                writer.writerow(row)


def main():
    if len(sys.argv) < 3:
        print("Usage: python additional_cleaning.py <input csv file> <output csv file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        clean_data(csv_file, output_file)
        print("Data cleaning complete. Output written to", output_file)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(2)


if __name__ == "__main__":
    main()
