# Amazon Web Services :cloud:

![AWS](https://img.shields.io/badge/AWS-%23FF9900.svg?style=for-the-badge&logo=amazon-aws&logoColor=white)

[AWS CLI Command Reference](https://docs.aws.amazon.com/cli/latest/index.html)

## EC2

[EC2 CLI Command Reference](https://docs.aws.amazon.com/cli/latest/reference/ec2/index.html)

Some useful tools for connecting to and interacting with ec2 instances;

- `ssh` for connecting to the instance
- `scp` and `rsync` for copying files to and from the instance
- `tmux` for persistent sessions
- `netrw` to remotely edit files and browse directories in `vim`

### Snapshots

aws ec2 describe-snapshots --owner-ids self --query 'Snapshots[*].{ID:SnapshotId,VolumeID:VolumeId,StartTime:StartTime,State:State,Progress:Progress,VolumeSize:VolumeSize,Description:Description}' --output table

```sh
echo ./update-aws-ssh.sh "$(aws ec2 describe-instances --instance-ids "$instance_id" --query 'Reservations[0].Instances[0].PublicDnsName' --output text)"
```

## S3

[S3 CLI Command Reference](https://docs.aws.amazon.com/cli/latest/reference/s3/index.html)

- s3 uri format `s3://bucketname/uri`

> [!WARNING]
> Object is of storage class GLACIER. Unable to perform download operations on GLACIER objects. You must restore the object to be able to the perform operation.

If you see the above warning when trying to make recursive calls to aws objects that have already been restored, rerun the command using the `--force-glacier-transfer` flag. There is no real error, this is just the way the aws cli handles the situation.

```sh
aws s3 cp s3://bucketname/uri/ . --recursive --force-glacier-transfer
```

Check out this [issue](https://github.com/aws/aws-cli/issues/1699) for more info.

```sh
# copy files to/from s3 bucket
aws s3 sync <source> <destination>

# copy all files of a certain type to local
aws s3 cp s3://bucketname/uri/*.fastq.gz .
```

- [AWS CLI Command Reference](https://docs.aws.amazon.com/cli/latest/index.html)
- [AWS S3 CLI Command Reference](https://docs.aws.amazon.com/cli/latest/reference/s3/index.html)
- s3 uri format `s3://bucketname/uri`

```sh
# copy files to/from s3 bucket
aws s3 sync <source> <destination>

# copy all files of a certain type to local
aws s3 cp s3://bucketname/uri/*.fastq.gz .
```

## Other CLI Tools

### `ssh`

edit `~/.ssh/config`

- add the location of the private key file so you don't have to specify it with `-i` each time you connect
- add the `User` and `Hostname` fields so you don't have to specify them each time

```sh
# example ~/.ssh/config
Host aws
    Hostname ec2-3-235-0-0.compute-1.amazonaws.com
    User ubuntu
    IdentityFile ~/.ssh/aws.pem
```

### `scp`

- "secure copy" files between computers using `ssh`
- this means the config from `~/.ssh/config` is used

```sh
# copy file to home dir of remote aws instance
scp localfile.txt aws:~

# copy remote file to pwd
scp aws:~/remotefile.txt .
```

### `rsync`

- useful for keeping all files in a directory up-to-date

```sh
# copy all files in a directory to a remote server
rsync -avz --progress /path/to/local/dir/ aws:/path/to/remote/dir/
```

### `tmux`

Weird stuff can happen with "nested" sessions over `ssh`.
If you want to attach to a tmux session on a remote server,
you need to use the `-t` flag since `tmux` is not a login shell.

```sh
ssh aws             # works
ssh aws tmux a      # huh?
ssh aws -t tmux a   # ok
```

### `vim`

Once you have ssh configured, you can use vim to edit files remotely thanks to
the `netrw` plugin that comes shipped with `vim`.

````sh
vim scp://aws/remote/path/to/file
Copy current vim buffer to remote server

``` vim
:!scp % aws:~/path/to/remote/file
````

Inide vim, you can use the `:Explore` command to browse the remote server.
