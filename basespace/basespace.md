### BaseSpace CLI (`bs`)
- [BaseSpace CLI Command Reference](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview)
``` sh
# install on macOS
# https://github.com/basespace/homebrew-basespace
brew tap basespace/basespace
brew install bs-cli

# then authenticate with your BaseSpace credentials
bs authenticate
```
| command| description |
|---|---|
`accession`|Register a new item
`add`|Add lane or workflow settings
`archive`|Archive entities
`authenticate`|Make an authentication request to BaseSpace (aliases: auth)
`await`|Await completion or status change
`children`|Retrieve child entities (aliases: child)
`clear`|Clear lane or workflow settings
`content`|Show all files in an entity (aliases: contents, dir)
`copy`|Copy an entity to a new location (aliases: duplicate)
`create`|Create an entity
`delete`|Move entities to your trash (aliases: rm, move-to-trash)
`download`|Download files from an entity
`empty`|Empty a store
`export`|Export lane or workflow settings
`get`|Get detailed information about an entity (aliases: info, inspect)
`header`|List headers for an entity (aliases: headers)
`history`|Retrieve account activity records for a user or workgroup
`import`|Import lane or workflow settings
`kill`|Abort entities
`launch`|Execute a task
`link`|Get a direct link for an entity
`list`|List and filter entities (aliases: filter, list-all)
`load`|Load into your environment
`logs`|Retrieve log files (aliases: log)
`rename`|Rename entities
`restore`|Restore items
`revoke`|Invalidate a resource (aliases: expire)
`seqstats`|Sequencing stats
`set`|Set properties in an entity
`translate`|Translate v1 <-> v2 entity IDs
`unarchive`|Restore entities from archival storage
`unlock`|Unlock a locked entity
`update`|Update entities
`upload`|Upload files to an entity
`whoami`|Get information about selected user/configuration
`yield`|Return yield information for an entity (aliases: yields)

#### `bs list project`

|Name|Id|
|---|---|
| HEL_H3K27ac_ChIPSeq1 |381715669|
| Untitled from 230306_NS500289_1231_AHW5Y3BGXN |383461079|
| CBP_P300 |388989668|



``` sh
bs list biosample > biosamples.txt
bs list biosample --project-id 381715669
bs list biosample --project-id 383461079
bs list biosample --project-id 388989668
```

