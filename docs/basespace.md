<!-- markdownlint-disable MD013 -->

# BaseSpace and the `bs` command-line interface

Read the [documentation](https://developer.basespace.illumina.com/docs).

Check out [examples](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-examples).

## Setup

Install on Linux or macOS:

```sh
sudo wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /usr/local/bin/bs
```

> [!CAUTION]
> If you want to install the binary somewhere other than `/usr/local/bin`,
> change the path in the command.

or install on macOS using Homebrew:

```sh
brew tap basespace/basespace && brew install bs-cli
```

After installation, authenticate with your BaseSpace credentials:

```sh
bs authenticate
```

Finally, run 'bs whoami' to verify that you are authenticated:

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

## Examples

### Downloading files

`bs download run --id 281093846`

```sh
ubuntu@ip-172-31-32-180:~/kalay/run4$ ls -p

Config/                      RTARead1Complete.txt
Data/                        RTARead2Complete.txt
InstrumentAnalyticsLogs/     RTARead3Complete.txt
InterOp/                     Recipe/
KB_ChIPseq_4_281093846.json  RunCompletionStatus.xml
Logs/                        RunInfo.xml
RTAComplete.txt              RunParameters.xml
RTAConfiguration.xml         SampleSheet.csv
RTALogs/
```

## Data

NextSeq run file downloaded from basespace and all runs are demultiplexed.

BaseSpace

| run date | run name (basespace) | size (GB) | %Q30  | %PF   | instrument | flow cell ID |
| -------- | -------------------- | --------- | ----- | ----- | ---------- | ------------ |
| 20230215 | HEL_H3K27ac_ChIPSeq1 | 27        | 94.49 | 92.21 | NB551203   | HC32JBGXN    |
| 20230306 | HEL_STAT5_H3K27ac_2  | 37        | 92.74 | 89.86 | NS500289   | HW5Y3BGXN    |
| 20230519 | CBP_p300             | 44        | 90.06 | 87.69 | NS500289   | HHLJNBGXT    |
| 20240522 | KB_ChIPseq_4         | 10        | 94.47 | 93.50 | NB551203   | H7TMNAFX7    |
| 20240617 | KB_ChIPseq_5         | 10        | ?     | ?     | ?          | ?            |

```plaintext
ubuntu@ip-172-31-32-180:~/kalay$ bs list runs
+---------------------------------+-----------+----------------------+----------+
|              Name               |    Id     |    ExperimentName    |  Status  |
+---------------------------------+-----------+----------------------+----------+
| 230215_NB551203_0628_AHC32JBGXN | 253069839 | HEL_H3K27ac_ChIPSeq1 | Complete |
| 230306_NS500289_1231_AHW5Y3BGXN | 254559306 | HEL_STAT5_H3K27ac_2  | Complete |
| 230519_NS500289_1262_AHHLJNBGXT | 258955746 | CBP_p300             | Complete |
| 240522_NB551203_0777_AH7TMNAFX7 | 281093846 | KB_ChIPseq_4         | Complete |
| 240617_NB551203_0782_AH7NV5AFX7 | 281976710 | Kalay5               | Failed   |
+---------------------------------+-----------+----------------------+----------+
```

Download runs with `bs download run --id "$run_id"`

bs download run --id 253069839
bs download run --id 254559306
bs download run --id 258955746
bs download run --id 281093846
bs download run --id 281976710

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
