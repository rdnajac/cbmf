import toml

# Read the TOML file
with open('palomero.toml', 'r') as file:
    data = toml.load(file)

# Extract author name
author_name = data.get('author', {}).get('name', 'Unknown Author')

# Extract publications
publications = data.get('publications', [])

# Generate Markdown content
markdown_content = f"# Publications by {author_name}\n\n"
for publication in publications:
    title = publication.get('title', 'No Title')
    journal = publication.get('journal', 'No Journal')
    year = publication.get('year', 'No Year')
    url = publication.get('url', '#')
    
    markdown_content += f"## {title}\n"
    markdown_content += f"*Journal:* {journal}\n"
    markdown_content += f"*Year:* {year}\n"
    markdown_content += f"*Link:* [{url}]({url})\n\n"

# Write the Markdown content to a file
with open('publications.md', 'w') as file:
    file.write(markdown_content)

print("Markdown file generated: publications.md")

