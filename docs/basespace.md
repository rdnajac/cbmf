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

### Scratch

Get the size of a run:
bs get run --id 281976710 --TotalSize
bs list runs -F Name -F QcStatus -F TotalSize -F AppSession.Application.Name

## bcl2fastq

```sh
bcl2fastq --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --auto-set-to-zero-barcode-mismatches --find-adapters-with-sliding-window --adapter-stringency 0.9 --mask-short-adapter-reads 35 --minimum-trimmed-read-length 35 -R run5 --sample-sheet ./run5/SampleSheet.csv -o ./run5_fastq
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
