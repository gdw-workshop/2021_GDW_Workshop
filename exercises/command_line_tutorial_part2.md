## Command-line clinic, continued 

### Paths: absolute, relative, and the $PATH

This morning, Bob introduced you to the concept of a path in unix, and explained how paths relate to the heirarchical filesystem.  Let's briefly revisit that topicby discussing absolute paths, relative paths, and the $PATH.  

#### absolute paths

Absolute paths start at the root of the filesystem.  In other words, absolute paths start with `/`.  `/Users/gdw/Desktop` is an absolute path. 

**Absolute paths always refer to the same place, no matter your present working directory**

#### relative paths

Relative paths don't start with an `/`.   `./Desktop` and `Desktop` and `../Desktop`  are all relative paths.  What they refer to depends on where you are. 

**A relative path's meaning depends on your present working directory**

#### PATH

The PATH has a special meaning in unix environments.  The PATH is a list of directories where the shell will look for commands.  If you didn't have the PATH, you'd always have to type the full path (absolute or relative) of a command to run it, which would be annoying.

PATH is what's called an [environmental variable](https://en.wikipedia.org/wiki/Environment_variable).  Here are two ways to find your PATH:

- Run the `env` command in the terminal.  Do you see your $PATH?
- Run `echo $PATH` in your terminal.  

Can you recognize what other environmental variables mean?

#### The PATH and the `which` command

The which command will you tell you where in the filesystem a command is, or more specically whether a command exists in your PATH.  To see a description of this command, run:

```
man which
```

Not to be confused with [manwich](https://manwich.com/) ![manwich](./manwich.jpg)

The `man` (manual) command is a great way to get more information about built-in linux commands.  You can page through a manual page using the up and down arrows or the space bar.  Press q to exit.  

OK, back to `which`.  All of the commands you run, even basic ones like `ls` and `cd`, are just a file somewhere in your path.  To find out where they are, run:

```
which ls
which cd
```

The PATH is a convience, so that instead of always typing:

```
/bin/ls
```

You can just type `ls`.  

Other software that's not a core part of the operating system, for instance the blastn and bowtie2 aligners, can also be found by `which`:

```
which blastn
which bowtie2
```

Even bash itself is a command in your PATH:
```
which bash
```

Note that all these commands are in your PATH.   Let's try something.  Let's change your PATH:

```
cd
export PATH=/Users/gdw
```

What do you think the impact of this change will be?  Try this:

```
env
```

How about these two commands?

```
ls
/bin/ls
../../bin/ls
../bin/ls
```

What's going on?  Why do these commands work or not work?

### Installing software from the command line

Installing software in a linux command-line environment can be a roadblock to beginners.  Let's practice installing a bioinformatics tool from the command line.

We'll install jellyfish, a tool for counting [k-mers](https://en.wikipedia.org/wiki/K-mer).  The main webpage for Jellyfish can be found [here](https://github.com/gmarcais/Jellyfish).  If you read that page, you will find installation instructions that direct you to the [releases page](https://github.com/gmarcais/Jellyfish/releases).

Here you find there are two options:

1. You can download a so-called pre-compiled binary.  Pre-compiled binaries are ready to go, but you must download a file that is matched to your operating system.  Since we are working on mac laptops, you'd want to download the jellyfish-macosx file from github

2. You can download the source code for the software and compile it to create a binary program file yourself. Sometimes this is the only option, for example for [minimap2](https://github.com/lh3/minimap2/releases), a program we'll use later this week.  


Let's first try downloading the pre-compiled binary.  Download the file named jellyfish-macosx
``` 
cd
curl -OL https://github.com/gmarcais/Jellyfish/releases/download/v1.1.12/jellyfish-macosx
```

`curl` is a program that will download files from the internet.

Type `ls -lh` to confirm that you have downloaded a file to the pwd named jellyfish-macosx that has a size of 434 kb.  The `-h` option to `ls` outputs file sizes in a **h**uman readable format and the `-l` option outputs a more detailed "**l**ong" listing.

Now that we have the file, let's try to run it.  That was the point of downloading it, after all.  Try typing:

```
jellyfish-macosx
```

You should have gotten an error like this:

```
-bash: jellyfish-macosx: command not found
```

What do you think is going on?

<br><br><br><br> <br><br><br><br> <br><br><br><br> <br><br><br><br> 

One issue is that the file jellyfish-macosx is not in your $PATH.  Or more accurately, jellyfish-macosx _is_ in your home directory, and your home directory is not in your $PATH.  Let's confirm this by using the which command again:

```
which jellyfish-macosx
```

To fix this situation, we could either move the jellyfish file into a directory in your PATH, or could add your home directory to your $PATH.  The first option would manifest as:

```
# option 1: move jellyfish to a directory already in your $PATH
sudo mv jellyfish-macosx /usr/local/bin
```

Note that you had to use the sudo (**s**uper **u**ser **do**) command to move jellyfish to /usr/local/bin.  This is because you lacked the necessary permissions to move a file into that directory.  Running commands prepended by sudo is the same as doing something with Administrator priveleges (so be careful!).

The second option would be:
```
# option 2: add ~ to your $PATH
export PATH=$PATH:/Users/gdw
```

Here, you are changing the $PATH the way you are supposed to, by appending a new directory to the $PATH instead of overwriting it as we did above.  So bash interprets `PATH=$PATH:/Users/gdw` as "assign to the variable named PATH the current value of the variable PATH plus ":/Users/gdw" 

After either (or both) of these, jellyfish-macosx should now be in your $PATH.  let's see what happens:

```
which jellyfish-macosx
jellyfish-macosx
```

You should see an error like this:

```
-bash: /usr/local/bin/jellyfish-macosx: Permission denied
```

Arrgh!  This is why people get frustrated with the command line.  Let's power through!

#### Permissions

This gets us to another common pitfall for linux beginners: [file permissions](http://linuxcommand.org/lc3_lts0090.php).

Every file and every directory in linux has a set of permissions that tell whether the file or directory can be read, written, or executed.  

Your user lacked write permissions for `/usr/local/bin`, which is why you couldn't move a file into that directory, and had to use the `sudo` command, which gives you super privileges.  

Similarly, the jellyfish-macosx file that you downloaded did not arrive with executable permissions. And in order to be run as a command, a file needs executable permissions.  

Let's check the permissions on jellyfish:

```
ls -l /usr/local/bin/
```

When you list a file using the `-l` option, the first part of the output for each file is the permissions, which look like this:

```
-rw-r--r--
```

These permissions have this meaning: ![permissions](./file_permissions.png)

You can see that this file indeed lacks e**x**ecutable permissions.  We can change that by running:

```
chmod +x /usr/local/bin/jellyfish-macosx
```

Now try running:

```
jellyfish-macosx
```

Phrew.  

Jellyfish got a little mad because you didn't supply enough arguments, but at least it ran!  That's a good sign.  You may also note that it outputs usage information about how to actually run it, which most well-written software should do.  You can also see that it has a `--help` option, which is another common feature of good software (sometimes the option is `-h` or `-help`).  

#### installing from source code

If you recall, we also could have downloaded the source code for jellyfish and compiled it to make our own executable program file.  Let's go through that quickly, just so you can see what that process looks like. 

First, download and unpack the source code

```
curl -OL https://github.com/gmarcais/Jellyfish/releases/download/v1.1.12/jellyfish-1.1.12.tar.gz
tar xvf jellyfish-1.1.12.tar.gz
```

This should have created a new directory named jellyfish-1.1.12.  We'll change into it and compile the software according to the instructions on the main github page: 

```
cd jellyfish-1.1.12
./configure --prefix=$HOME
make -j 4
```

Note that `$HOME` in the above code refers to the value of the environmental variable HOME, which is the path of your home directory.  When bash runs the above code, it will replace `$HOME` with the value `/Users/gdw`

You will see a bunch of computer sciencey things print to the terminal as you run this command.  There may be a warning, but there should be no errors, which would indicate that compilation failed.  This compilation process created a binary executable file in the bin directory.  Look at it by running:

```
ls bin
```

This process created a file named jellyfish that is the same size as the jellyfish-macosx file that we downloaded previously.  Let's try running it.  Note that this directory is not in your PATH, so we'll have to refer to it explicitly by running:

```
bin/jellyfish
```

You should see the same output as you saw previously. 

### Operating more efficiently in the command line environment

There are a number of simple tricks that will enable you to operate more efficiently in the command line environment.


#### Tab completion

[Tab completion](https://en.wikipedia.org/wiki/Command-line_completion) will allow you to type less by automatically filling in the names of files and commands when possible.  It works in 2 ways:

1. If you are at the beginning of a command line, it will complete using the already typed characters and the names of the commands available in your PATH.  

2. If you are in the middle of a command line, it will complete using the already typed characters and the names of files and directories in your pwd.  

If you hit `tab` one time, and there is an unambiguous matching command or file the could be completed, it will complete it as far as it can until there is ambiguity again.   If you hit `tab` twice, it will show you all the possibly matching files or commands.

For instance, say you want to run bowtie2.  You could type out bowtie2, or you could:

##### command completion 

- enter `bowti` on the command line and hit `tab` once.  It should complete to `bowtie2`
- try hitting `tab` twice more at this point.  You will see the other commands that begin with bowtie2.
- enter `b` on the command line and hit `tab` twice.  You will see all the commands that begin with `b`.

##### command completion 

`cd` to get to your home directory.  Then type:

- `ls D` then hit `tab` twice.  You should see Desktop, Downloads, and Documents, which are all the directories that begin w/ D
- now type `ls De` and hit `tab` once.  This should complete to Desktop, since De is enough to distinguish Desktop from Downloads and Documents.

Becoming comfortable with tab completion will make your life much easier when typing commmands.

### Wildcards

Wildcards are somewhat related to tab completion in that they allow specified ambiguity.  When bash is interpreting a wildcard character, it will replace the word containing a wildcard with all the possible matching files or directories.  Some useful wildcards are:

- `*`     match any number of any character
- `?`     match one of any character
- `[123]` match one of the characters `1`, `2`, or `3` 

There are a number of wildcards that we won't go into, but you can read more [here](http://tldp.org/LDP/GNU-Linux-Tools-Summary/html/x11655.htm)

For example:

```
ls /usr/local/bin/b*
ls /usr/local/bin/*seq*
ls /usr/local/bin/bowtie2-align-?
ls /usr/local/bin/bowtie2-align-[sl]
```

We will use wildcards more during the workshop.

### Aliases

Creating command aliases can be a convenient shortcut.  They allow you to define a new command or to redine the meaning of a command, for instance:

```
alias 'ls=ls -G'
ls
alias 'll=ls -l'
ll
```

Want to know what the `-G` and `-l` options do for `ls`?  See `man ls`.

Now, with these aliases, when you type `ls`, the shell will replace `ls` with `ls -G`.  And when you type `ll`, it will replace `ll` with `ls -l`, which will in turn be expanded to `ls -G -l`


### Connecting to remote servers

#### ssh remote shell

We already used `ssh` to connect to a remote server when we were checking the resources available on the cctsi-104 server.  `ssh` provides a shell (like bash) on a remote server.   Let's try connecting again to the cctsi-104 server:

```
ssh gdw@cctsi-104.cvmbs.colostate.edu
```

Once logged in, use `ls` to see what files are in the gdw user's home directory on cctsi-104.  

#### sftp file transfer

`ssh` gives you a remote shell.  `sftp` is a command line tool that allows you to copy files back and forth from remote servers.  (`sftp` is the secure (encrypted) version of `ftp`, **f**ile **t**ransfer **p**rotocol).

Let's get those dataset (.fastq files) from the cctsi-104 server.  

```
sftp gdw@cctsi-104.cvmbs.colostate.edu
```

`sftp` has some limited shell-like functionality.  For instance, you can `cd`, `ls`, and `pwd` from within your sftp session.  By default, you will be logged into your home directory.  Let's double check that those fastq files are there and then use the `get` command within sftp to transfer them from the remote server to your local computer:

```
# within sftp sessions
pwd
ls 
# wildcards work within sftp
get *.fastq 
```

Now exit from the sftp session by typing `exit`

You should see that you transferred 3 fastq format sequence files to whatever directory you were in when you ran `sftp`.

`put` is the opposite of `get` in sftp.  You can use it to transfer files _to_ a remote server from the command line.


### The shell history, pipes, and `grep`

#### history

You might want to run a command again or just remember how you ran a command previously.  bash keeps track of the commands you ran and you can view a record of your previous commands using the `history` command.  Try running `history`.


#### pipes

One really useful feature of bash and similar shells is the ability to use the output of one command as the input of another command.  This process is called piping.  The symbol `|` in bash is called a pipe.  Here's a simple example of how you could pipe 2 commands together:

```
# show me the files in /usr/local/bin that begin with bowtie
ls /usr/local/bin/bowtie*
# use the wc -l command to count those files
ls /usr/local/bin/bowtie* | wc -l
```

#### grep

`grep` is a command that will search for a pattern in its input and report lines that match that pattern.  The input to `grep` could be a file or it could be input from a pipe, as in this example:

```
# show me the files in /usr/local/bin
ls /usr/local/bin
# in the output of that command, search for files that match the pattern "bowtie"
ls /usr/local/bin | grep bowtie
# you can pipe multiple commands together to create a "pipeline"
ls /usr/local/bin | grep bowtie | wc -l
```

Note that this last command did the same thing as `ls /usr/local/bin/bowtie* | wc -l`.  In bash there are often many ways to do the same thing.



### bash scripting

Any file that you make executable that contains bash commands is a bash script. 

Let's make some simple scripts to demystify the process.  First, let's make a script that contains the command we use to login to the cctsi-104 remote server.  Open BBEdit, and add the following one line to a file: 

```
ssh gdw@cctsi-104.cvmbs.colostate.edu
```

Save the file to the Desktop, naming it ssh_104.  

Now find the file:

```
cd ~/Desktop
ls -l ssh_104
```

What would happen if you tried to run this script?  What command do we need to run to fix it?

<br><br><br><br> <br><br><br><br> <br><br><br><br> <br><br><br><br> 

That's right, it doesn't have executable permissions.  Let's add them using `chmod`

```
chmod +x ssh_104
```

Now we can run this script and it will execute whatever commands are listed in it.  In this case, it will run `ssh gdw@cctsi-104.cvmbs.colostate.edu`

Since ~/Desktop is not in our PATH, we have to specify the full path to the script using ./ (assumes pwd is ~/Desktop)

```
./ssh_104
```

Bash is also a full-fledged scripting language, and your script can include flow control (if/else statements, loops), and variables.  Let's make a script that includes a for loop and a variable.  Open BBEdit and enter this into a file:

```
# a for loop 
for bowtie_file in /usr/local/bin/bowtie*
do

	# $bowtie_file is a variable whose value changes each loop.
   echo "I found this bowtie file in /usr/local/bin: $bowtie_file"
done
```

Save this file to the Desktop, naming it bowtie_script or something like that, give it executable permissions, and run it.  Did it work?


