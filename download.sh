#!/bin/bash

wget -O test_data.tar.gz "https://www.dropbox.com/scl/fi/9dmjdg6f6kqpkgbxm87t7/test_data.tar.gz?rlkey=qdvk54r3dq7ib5t46pnuvjeka&st=b7fc5m5c&dl=0"
wget -O 1000G_hg38.zip "https://www.dropbox.com/scl/fi/d0xm2ht5scduspucmr2qb/1000G_hg38.zip?rlkey=niv01mpleiqjosvcljgguf6t9&st=9lxnh2tf&dl=0"
tar -xzvf test_data.tar.gz
unzip 1000G_hg38.zip
