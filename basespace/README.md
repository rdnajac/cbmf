## BaseSpace CLI: `bs`

Basespace CLI is a command-line interface for Illumina BaseSpace. It allows you to interact with BaseSpace from the command line, enabling you to automate workflows and integrate BaseSpace with other tools.

## Installation


Ubuntu:

```sh
sudo wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /usr/local/bin/bs


```

### macOS

To install the BaseSpace CLI on macOS, use Homebrew:

```sh
brew tap basespace/basespace
brew install bs-cli
```

### Authenticate

After installation, authenticate with your BaseSpace credentials:

```sh
bs authenticate
```

### Command Reference

The following table provides a quick reference for the most commonly used BaseSpace CLI commands.

| Command        | Description                                                       |
| -------------- | ----------------------------------------------------------------- |
| `accession`    | Register a new item                                               |
| `add`          | Add lane or workflow settings                                     |
| `archive`      | Archive entities                                                  |
| `authenticate` | Make an authentication request to BaseSpace (aliases: auth)       |
| `await`        | Await completion or status change                                 |
| `children`     | Retrieve child entities (aliases: child)                          |
| `clear`        | Clear lane or workflow settings                                   |
| `content`      | Show all files in an entity (aliases: contents, dir)              |
| `copy`         | Copy an entity to a new location (aliases: duplicate)             |
| `create`       | Create an entity                                                  |
| `delete`       | Move entities to your trash (aliases: rm, move-to-trash)          |
| `download`     | Download files from an entity                                     |
| `empty`        | Empty a store                                                     |
| `export`       | Export lane or workflow settings                                  |
| `get`          | Get detailed information about an entity (aliases: info, inspect) |
| `header`       | List headers for an entity (aliases: headers)                     |
| `history`      | Retrieve account activity records for a user or workgroup         |
| `import`       | Import lane or workflow settings                                  |
| `kill`         | Abort entities                                                    |
| `launch`       | Execute a task                                                    |
| `link`         | Get a direct link for an entity                                   |
| `list`         | List and filter entities (aliases: filter, list-all)              |
| `load`         | Load into your environment                                        |
| `logs`         | Retrieve log files (aliases: log)                                 |
| `rename`       | Rename entities                                                   |
| `restore`      | Restore items                                                     |
| `revoke`       | Invalidate a resource (aliases: expire)                           |
| `seqstats`     | Sequencing stats                                                  |
| `set`          | Set properties in an entity                                       |
| `translate`    | Translate v1 <-> v2 entity IDs                                    |
| `unarchive`    | Restore entities from archival storage                            |
| `unlock`       | Unlock a locked entity                                            |
| `update`       | Update entities                                                   |
| `upload`       | Upload files to an entity                                         |
| `whoami`       | Get information about selected user/configuration                 |
| `yield`        | Return yield information for an entity (aliases: yields)          |

## Examples

### `bs list runs`

Running `bs list runs` will list all runs in your BaseSpace account.

```plaintext
+---------------------------------+-----------+----------------------+----------+
|              Name               |    Id     |    ExperimentName    |  Status  |
+---------------------------------+-----------+----------------------+----------+
| 230215_NB551203_0628_AHC32JBGXN | 253069839 | HEL_H3K27ac_ChIPSeq1 | Complete |
| 230306_NS500289_1231_AHW5Y3BGXN | 254559306 | HEL_STAT5_H3K27ac_2  | Complete |
| 230519_NS500289_1262_AHHLJNBGXT | 258955746 | CBP_p300             | Complete |
| 240522_NB551203_0777_AH7TMNAFX7 | 281093846 | KB_ChIPseq_4         | Complete |
+---------------------------------+-----------+----------------------+----------+
```
## Download the bcl files for a run

```sh
bs download run --id 281093846
```
```
