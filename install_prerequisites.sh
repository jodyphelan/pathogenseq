BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir bin
# bcftools
wget https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2 \
&& tar -xvf bcftools-1.6.tar.bz2 \
&& cd bcftools-1.6 \
&& ./configure \
&& make \
&& cp bcftools ../bin/ \
&& cd htslib-1.6 \
&& make  \
&& cp bgzip tabix ../../bin/ \
&& cd ../../

cd $BASE_DIR
# samtools
wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 \
&& tar -xvf samtools-1.7.tar.bz2 \
&& cd samtools-1.7 \
&& ./configure \
&& make \
&& cp samtools ../bin/ \
&& cd ../

cd $BASE_DIR
# htsbox
git clone https://github.com/lh3/htsbox.git \
&& cd htsbox  \
&& make \
&& cp htsbox ../bin/ \
&& cd ../

cd $BASE_DIR
# bwa
git clone https://github.com/lh3/bwa.git \
&& cd bwa \
&& make \
&& cp bwa ../bin/ \
&& cd ../

cd $BASE_DIR
# minimap2
git clone https://github.com/lh3/minimap2.git \
&& cd minimap2 \
&& make \
&& cp minimap2 ../bin/ \
&& cd ../

cd $BASE_DIR
# examl
git clone https://github.com/stamatak/ExaML.git \
&& cd ExaML/parser \
&& make -f Makefile.SSE3.gcc \
&& cp parse-examl ../../bin/ \
&& cd ../examl \
&& make -f Makefile.OMP.AVX.gcc && cp examl-OMP-AVX ../../bin/examl || (make -f Makefile.OMP.SSE3.gcc && cp examl-OMP ../../bin/examl) \
&& cd ../../

cd $BASE_DIR
# raxml
git clone https://github.com/stamatak/standard-RAxML.git \
&& cd standard-RAxML/ \
&& (make -f Makefile.AVX2.PTHREADS.gcc && cp raxmlHPC-PTHREADS-AVX2 ../bin/raxmlHPC) || (make -f Makefile.AVX2.PTHREADS.gcc && cp raxmlHPC-PTHREADS-AVX ../bin/raxmlHPC) || (make -f Makefile.SSE3.PTHREADS.gcc && cp raxmlHPC-PTHREADS-SSE3 ../bin/raxmlHPC) \
&& cd ../

cd $BASE_DIR
pip install -r requirements.txt

echo "echo \"export PATH=$PWD:$PATH\" >> ~/.bashrc"
