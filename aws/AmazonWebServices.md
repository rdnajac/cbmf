# Amazon Web Services

## EC2

## S3

### `aws`

- [AWS CLI Command Reference](https://docs.aws.amazon.com/cli/latest/index.html)
- [AWS S3 CLI Command Reference](https://docs.aws.amazon.com/cli/latest/reference/s3/index.html)
- s3 uri format `s3://bucketname/uri`

``` sh
# copy files to/from s3 bucket
aws s3 sync <source> <destination>

# copy all files of a certain type to local
aws s3 cp s3://bucketname/uri/*.fastq.gz .
```

## Remote access

- Connect to and interact the ec2 instance using `ssh` and `scp`
- Bonus points if you `netrw` to edit files remotely with `vim`
- TODO `tmux` for persistent sessions
- TODO `rsync` for keeping directories in sync
- TODO add section on ftp and downloading from Azenta

### `aws`

- [AWS CLI Command Reference](https://docs.aws.amazon.com/cli/latest/index.html)
- [AWS S3 CLI Command Reference](https://docs.aws.amazon.com/cli/latest/reference/s3/index.html)
- s3 uri format `s3://bucketname/uri`

``` sh
# copy files to/from s3 bucket
aws s3 sync <source> <destination>

# copy all files of a certain type to local
aws s3 cp s3://bucketname/uri/*.fastq.gz .
```

### `ssh`

edit `~/.ssh/config`
- add the location of the private key file so you don't have to specify it with `-i` each time you connect
- add the `User` and `Hostname` fields so you don't have to specify them each time
``` sh
# example ~/.ssh/config
Host aws
    Hostname ec2-3-235-0-0.compute-1.amazonaws.com
    User ubuntu
    IdentityFile ~/.ssh/aws.pem
```

### `scp`

- "secure copy" files between computers using `ssh`
- this means the config from `~/.ssh/config` is used
``` sh
# copy file to home dir of remote aws instance
scp localfile.txt aws:~

# copy remote file to pwd
scp aws:~/remotefile.txt .
```

### `rsync`

- useful for keeping all files in a directory up-to-date

``` sh
# copy all files in a directory to a remote server
rsync -avz --progress /path/to/local/dir/ aws:/path/to/remote/dir/
```

### `tmux`

Weird stuff can happen with "nested" sessions over `ssh`. If you want to attach to a tmux session on a remote server, you need to use the `-t` flag since `tmux` is not a login shell.
``` sh
ssh aws             # works
ssh aws tmux a      # huh?
ssh aws -t tmux a   # ok
```

#### `vim`

Once you have ssh configured, you can use vim to edit files remotely thanks to the `netrw` plugin that comes shipped with `vim`.
``` sh
vim scp://aws/remote/path/to/file
Copy current vim buffer to remote server

``` vim
:!scp % aws:~/path/to/remote/file
```
