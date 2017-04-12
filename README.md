# BYU ACME Labs #

The file `README.md` produces a formatted introduction to this repository, displayed on the **Overview** and **Source** pages.
Feel free to edit this file and customize the repository.
Read [Bitbucket's Markdown tutorial](https://bitbucket.org/tutorials/markdowndemo) to learn _Markdown_, the language that the file is written in.

## Repository Organization ##

Labs are either submitted as `.py` files (Python) or `.ipynb` files (Jupyter Notebook.
These files should be placed in the main directory of the repository.

Some labs are graded with an automated test driver.
The TA will place feedback from the test drivers in the `Feedback/` folder.
Apart from these feedback files, files with an extension other than `.py` or `.ipynb` are not tracked by this online repository.

## Git ##

This website is a _git repository_, an online storage space for code and other small files.
_Git_ is the underlying software that manages updates between this online repository and the copies (or _clones_) of the repository stored locally on computers.
Most computers already come with git installed, but it can be downloaded [here](http://git-scm.com/downloads) if needed.

#### Cloning ####

To get started, open up a shell application.
This is usually called _Terminal_ on Linux or Mac, or _GitBash_ on Windows (GitBash comes with the Windows git installation).
Navigate to the folder where you want to store the files from this repository.
Use **cd <directory>** to move to the specified <directory> and **pwd** to print the working directory.
Then run **git clone <repo_url>**.
You can find <repo_url> by clicking **Clone** on the left, under **ACTIONS**.
This action creates a new folder that serves as a copy of this repository.

#### Using Git ####

Git manages the history of a file system through a series of _commits_, or checkpoints.
Use **git status** to see the files that have been changed since the last commit.
These changes are then moved to the _staging area_, a list of files to save during the next commit, with **git add <filename(s)>**.
Save the changes in the staging area with **git commit -m "A brief message describing the changes"**.

All of these commands are done within a clone of the repository, stored somewhere on a computer.
This repository must be manually synchronized with the online repository via two other git commands: **git pull origin master**, to pull updates from the web to the computer; and **git push origin master**, to push updates from the computer to the web.

#### Common Git Commands ####

| Command                     | Explanation                                    |
|:------------------------------------------------|:---------------------------|
| git status                  | Display the staging area and untracked changes.|
| git pull origin master      | Pull down changes from the online repository.  |
| git push origin master      | Push up changes to the online repository.      |
| git add <filename(s)>       | Add a file or files to the staging area.       |
| git add -u                  | Add all of the modified, previously tracked files to the staging area.|
| git commit -m "<message>"   | Save the changes in the staging area with a given message.|
| git checkout -- <filename>  | Revoke the changes made since the last commit on a file that is not in the staging area. |
| git reset -- <filename>     | Remove a file from the staging area.           |
| git diff <filename>         | See the changes made on an unstaged file since the last commit.|
| git diff --cached <filename> | See the changes made on a staged file since the last commit.|

#### Example Work Session ####

Short version:
```
#!bash
$ cd Desktop/Vol1A/
$ git pull origin master                           # Pull updates.
$ touch python_intro.py                            # Make changes.
$ git add python_intro.py                          # Track changes.
$ git commit -m "Made a file for the first lab."   # Commit changes.
$ git push origin master                           # Push updates.
```

Long version:
```
#!bash
# Navigate to the clone of the repository.
$ cd Desktop/Vol1A

# Pull any updates from the online repository (such as TA feedback).
$ git pull origin master
From https://bitbucket.org/byuacmeta/template
 * branch            master     -> FETCH_HEAD
Already up-to-date.

# Work on the labs. For example, create a new file for the first lab.
$ touch python_intro.py
$ git status
On branch master
Your branch is up-to-date with 'origin/master'.
Untracked files:
  (use "git add <file>..." to include in what will be committed)

	python_intro.py

nothing added to commit but untracked files present (use "git add" to track)

$ git add python_intro.py 
$ git status
On branch master
Your branch is up-to-date with 'origin/master'.
Changes to be committed:
  (use "git reset HEAD <file>..." to unstage)

	new file:   python_intro.py

# Commit the changes to the repository with an informative message.
$ git commit -m "Made a file for the first lab"
[master fed9b34] Made a file for the first lab
 1 file changed, 0 insertions(+), 0 deletions(-)
 create mode 100644 python_intro.py

# Push the changes to the online repository.
$ git push origin master
Counting objects: 3, done.
Delta compression using up to 2 threads.
Compressing objects: 100% (2/2), done.
Writing objects: 100% (3/3), 327 bytes | 0 bytes/s, done.
Total 3 (delta 0), reused 0 (delta 0)
To https://byuacmeta@bitbucket.org/byuacmeta/template.git
   5742a1b..fed9b34  master -> master

# The changes have been saved and the online repository updated.
$ git status
On branch master
Your branch is up-to-date with 'origin/master'.
nothing to commit, working directory clean
```# Math_Analysis_16
