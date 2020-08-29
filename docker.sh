# the gitlab testing environment relies on docker image.
# If the test failed online and the error is hard to debug, it
# might be useful to setup a local docker container to run the 
# test within in order to reproduce the error. In this case, follow
# the instruction below.


# pull docker image with python 3.7 (debian stretch)
docker pull python:3.7-stretch

# run bash shell on the docker image
docker run -it python:3.7-stretch bash

# this will create a docker image in which you will need to run
# the commands in .gitlab-ci.yml in order to install necessary
# packages, then you can run the test `python test.py`


