# Easy-Web-BLAST

## Introduction

Easy-Web-BLAST (EWB) is a single CGI script which can be deployed by a user to make a custom BLAST web page on a Unix server. It's key features are:

* one simple configuration file 
* automatically detects all the installed BLAST databases
* automagically chooses the correct BLAST tool to use (blastp, tblastn, ...)
* uses the latest BLAST+ toolkit


## Installation

There are a few different ways in which EWB can be installed. Some choices require root access to your server, and some can be done purely at the user level.

### Apache2

    cd $HOME/public_html
    tar zxvf easy-web-blast.tar.gz
    mv easy-web-blast myname


Then go to: http://example.com/~user/myname/

### Python Server

FIXME

## Configuration

You must edit the blast.conf file to set up the BLAST server:
* exe_path
  * Folder containing the BLAST+ binaries
* db_path
  * Folder containing the formatted BLAST indices
* admin_name
  * Name of person responsible for your web site
* admin_email
  * Their email address
* name
  * Title of the web page eg. "John's Parasite BLAST page"
* logo
  * An image to display eg. the institute logo. Can be a file in the EWB folder, or a global URL.

## FAQ

* How do I restrict which databases in the folder are visible? 
  * You can use symlinks to make a virtual database folder. First make a folder in the EWB folder called "db" and inside DB make symlinks to the indices on your server. Remember to change the db_path variable in blast.conf to be "./db" now instead.

## Authors

* Torsten Seemann <torsten.seemann@monash.edu>
* David Powell <david.powell@monash.edu>

## Reference

https://github.com/Victorian-Bioinformatics-Consortium/easy-web-blast
