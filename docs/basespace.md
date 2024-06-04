# BaseSpace CLI: `bs`

Read the [documentation](https://developer.basespace.illumina.com/docs).

## Setup

Skip ahead to [examples](#examples) if you already `bs` installed and authenticated.

### Installation

Ubuntu:

```sh
sudo wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /usr/local/bin/bs
```

> [!CAUTION]
> This command will download the BaseSpace CLI executable to `/usr/local/bin/bs`.
> If you want to install elsewhere, change the path in the command.

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

### Verify

Run 'bs whoami' to verify that you are authenticated:

```plaintext
+----------------+----------------------------------------------------+
| Name           | Ryan Najac                                         |
| Id             | ########                                           |
| Email          | rdn2108@cumc.columbia.edu                          |
| DateCreated    | 2021-07-13 15:29:51 +0000 UTC                      |
| DateLastActive | 2024-06-03 18:59:47 +0000 UTC                      |
| Host           | https://api.basespace.illumina.com                 |
| Scopes         | READ GLOBAL, CREATE GLOBAL, BROWSE GLOBAL,         |
|                | CREATE PROJECTS, CREATE RUNS, START APPLICATIONS,  |
|                | MOVETOTRASH GLOBAL, WRITE GLOBAL                   |
+----------------+----------------------------------------------------+
```

ubuntu@ip-172-31-32-180:~$

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

## Download the bcl files for all runs

```sh
# capture the output following the second 
for run in $(bs list runs --json | jq -r '.[].Id'); do
  bs download run $run --output-dir $run
done
one liner:
for run in $(bs list runs --json | jq -r '.[].Id'); do bs download run $run --output-dir $run; done
```

```

```

## Command Reference

The following table provides a quick reference for `bs` commands:

| Command        | Description                                                 |
| -------------- | ----------------------------------------------------------- |
| `accession`    | Register a new item                                         |
| `add`          | Add lane or workflow settings                               |
| `archive`      | Archive entities                                            |
| `authenticate` | Make an authentication request to BaseSpace (aliases: auth) |
| `await`        | Await completion or status change                           |
| `children`     | Retrieve child entities (aliases: child)                    |
| `clear`        | Clear lane or workflow settings                             |
| `content`      | Show all files in an entity (aliases: contents, dir)        |
| `copy`         | Copy an entity to a new location (aliases: duplicate)       |
| `create`       | Create an entity                                            |
| `delete`       | Move entities to your trash (aliases: rm, move-to-trash)    |
| `download`     | Download files from an entity                               |
| `empty`        | Empty a store                                               |
| `export`       | Export lane or workflow settings                            |
| `get`          | Get information about an entity (aliases: info, inspect)    |
| `header`       | List headers for an entity (aliases: headers)               |
| `history`      | Retrieve account activity records for a user or workgroup   |
| `import`       | Import lane or workflow settings                            |
| `kill`         | Abort entities                                              |
| `launch`       | Execute a task                                              |
| `link`         | Get a direct link for an entity                             |
| `list`         | List and filter entities (aliases: filter, list-all)        |
| `load`         | Load into your environment                                  |
| `logs`         | Retrieve log files (aliases: log)                           |
| `rename`       | Rename entities                                             |
| `restore`      | Restore items                                               |
| `revoke`       | Invalidate a resource (aliases: expire)                     |
| `seqstats`     | Sequencing stats                                            |
| `set`          | Set properties in an entity                                 |
| `translate`    | Translate v1 <-> v2 entity IDs                              |
| `unarchive`    | Restore entities from archival storage                      |
| `unlock`       | Unlock a locked entity                                      |
| `update`       | Update entities                                             |
| `upload`       | Upload files to an entity                                   |
| `whoami`       | Get information about selected user/configuration           |
| `yield`        | Return yield information for an entity (aliases: yields)    |
