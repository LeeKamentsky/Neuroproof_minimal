wget https://s3-us-west-2.amazonaws.com/sct-cloudmachines/new+folder/cilkplus-install.tar.gz
tar xvf cilkplus-install.tar.gz
ln -s $PWD/cilkplus-install/lib64/libcilkrts.so.0 $CONDA_PREFIX/lib/libcilkrts.so.0 
