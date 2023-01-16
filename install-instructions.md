# Clone from repository
git clone git@github.com:Mengzhendada/sidis.git # I forked to my git, but duane might update it

cd sidis 

#There are some git submodule things, which is difficult to clone proper path, so I just cloned the origin position

cd external

#There are three directories, bubble, cubature-cpp, mstwpdf, I just downloaded from Duane origin repository
git clone git@github.com:duanebyer/bubble.git
git clone git@github.com:duanebyer/cubature-cpp.git
git clone git@github.com:duanebyer/mstwpdf.git
#These three directories are supposed to be linked by git submodule

#Build using CMake
mkdir build
cd build
#make sure to use the correct apps, gcc/9.3.0, cmake/3.23.2 and root
source /group/solid/apps/root/root_v6.24.06.Linux-centos7-x86_64-gcc4.8/bin/thisroot.sh

#For me, even though gcc is updated to the latest version, g++ also latest, but not cc, so I have to specify which compiler to use in cmake

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/apps/gcc/9.3.0/bin/gcc -DCMAKE_CXX_COMPILER=/apps/gcc/9.3.0/bin/g++ ..

make
#make VERBOSE=1 for details

Now the executable file is in build/app/sidisgen

