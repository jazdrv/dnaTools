========================================================
some quick notes about the folder + setup scheme so far
========================================================

---------------------------------------------------
Directory Structure:
---------------------------------------------------

./LICENSE - open GL 3 license (got this from github)

./README.md - originally pushed this in using bitbucket, it's the project default
readme at GitHub

./bin - this is the place I am setting up for where the executables are to be found

./data - csv's and other types of data dumps

./docs - pdf's, etc

./env - workspace area (sort of what I see redux.bash was doing in that prep work
for zips/unzips/txt files and such) ... a separate space so we can separate it
from the real code directory; i've been using sym-links for some of these
things to point to stuff I have in data directory, etc

./sql - just sql files

./data - csv's, zips, etc 

./versions - old versions

./tests - place for unit tests (back burner for now)

./bin/redux - place for where the redux.bash (and supporting files) is

./bin/redux2 - place for the redux2.py script (and supporting files) is; i
haven't been editing redux2.py ... because I don't want to break anything. I've
been working with redux2z.py; my version needs db.py, lib.py; i have clades.py
in here and unpack.py ... but I haven't really studied that code yet

note: for my directory structure, i am doing a sym links like this in my bin folder:

bin jazdrv$ ls -l|grep ">"
lrwxr-xr-x   1 jazdrv  staff     20 Jan 12 09:01 clades2.py -> ./redux2/clades2z.py
lrwxr-xr-x   1 jazdrv  staff     16 Jan 11 11:11 redux.sh -> ./redux/redux.sh
lrwxr-xr-x   1 jazdrv  staff     19 Jan 11 12:54 redux2.py -> ./redux2/redux2z.py

note: I'm not tied to any one particular folder structure over another.
so feel free to change any of this.

---------------------------------------------------
Setup
---------------------------------------------------

I've been operating under the idea that a yaml config can be used. See:
./bin/redux2/config.yaml

And also one or two things that can be set as environment variables:

export REDUX_CONF="/Users/jazdrv/_prj/dnatools/bin/redux2/config.yaml"
export PYTHONPATH="$PYTHONPATH:/Users/jazdrv/_prj/dnatools/bin/redux2"
export REDUX_ENV="/Users/jazdrv/_prj/dnatools/env"
export REDUX_SQL="/Users/jazdrv/_prj/dnatools/sql"
export REDUX_PATH="/Users/jazdrv/_prj/dnatools/bin/redux2"
export REDUX_DATA="/Users/jazdrv/_prj/dnatools/data"

I don't think all these environment var's need setting like this though. They can be done
in the config file too in a better way. 

Again, I'm not tied to anything in the setup either. I can accomodate
multiple ideas -- whatever makes sense.



