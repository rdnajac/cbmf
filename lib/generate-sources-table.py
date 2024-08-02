import toml


# Read TOML file
def read_toml(file_path):
    with open(file_path, "r") as file:
        data = toml.load(file)
    return data


# Generate Markdown table
def generate_markdown(data):
    markdown = "# Sources\n\n"
    markdown += "Source list and software installation instructions.\n\n"
    markdown += "## Software\n\n"
    markdown += "| Software | Repository | Publication | DOI | Version |\n"
    markdown += "| --- | --- | --- | --- | --- |\n"

    footnotes = []
    footnote_counter = 1

    for software in data["software"]:
        name = software.get("name", "")
        repository = software.get("repository", "N/A")
        publication_name = software.get("publication_name", "N/A")
        doi = software.get("doi", "N/A")
        version = software.get("version", "N/A")

        if doi != "N/A":
            footnotes.append(f"[^{footnote_counter}]: {publication_name} ({doi})")
            publication_ref = f"[{publication_name}][^{footnote_counter}]"
            footnote_counter += 1
        else:
            publication_ref = publication_name

        markdown += (
            f"| {name} | {repository} | {publication_ref} | {doi} | {version} |\n"
        )

    markdown += "\n"
    markdown += "\n".join(footnotes)

    return markdown


# Save Markdown to file
def save_markdown(content, file_path):
    with open(file_path, "w") as file:
        file.write(content)


# Read TOML data
data = read_toml("sources.toml")

# Generate Markdown content
markdown_content = generate_markdown(data)

# Save Markdown content to README.md
save_markdown(markdown_content, "README.md")
