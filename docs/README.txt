========================================================
some quick notes about the folder + setup scheme so far
========================================================

---------------------------------------------------
Directory Structure:
---------------------------------------------------

./LICENSE - open GPL 3 license (got this from github)

./README.md - the project default readme at GitHub

./attic - some old things that might still be useful for reference

./bin - executables for current version

./data - csv's and other types of data dumps
./data/HaplogroupR - mirror of Hap-R data warehouse directory for zip files

./docs - pdf's, READMEs, instructions, etc

./examples - code snippets, reference examples

./sql - just sql files

./src - source code, where the programs live

./tests - place for unit tests (back burner for now)

./versions - old versions


---------------------------------------------------
Setup
---------------------------------------------------

1. Edit the config file. Really. Look at all of the variables near the top.

./src/config.yaml

2. Set the environment variable

. ./redux.env

3. Run the driver, which creates and populates the database basics. This is
an important step in the setup if it has never been run before

cd ${REDUX_PATH}
./redux.py -c

4. Download some data files - supply proper user credentials for
ftp. It's not important which set of or how many files you download.

sqlite3 variant.db > /tmp/f.out <<EOF
select fileNm from dataset where fileNm like '%/b38/%';
EOF

cd ../data/HaplogroupR
IFS_backup=$IFS
IFS=$(echo -en "\n\b")
for f in `head -10 /tmp/f.out`; do wget -nc -x -nH --user <user>@it2kane.org --password=<password> ftp://it2kane.org/"$f"; done
IFS=${IFS_backup}

5. Run the driver in test mode

./redux.py -t

<stay tuned for more progress>
