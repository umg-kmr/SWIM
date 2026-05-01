FROM python:3.14-slim

ENV DEBIAN_FRONTEND=noninteractive

EXPOSE 8888

# -----------------------------
# System dependencies
# -----------------------------
RUN apt-get update && apt-get install -y \
    build-essential htop libopenblas0 \
    gcc g++ gfortran make cmake bash-completion \
    libopenmpi-dev openmpi-bin \
    git curl wget ca-certificates vim nano sudo \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------
# Create non-root user
# -----------------------------
RUN useradd -m -s /bin/bash swim && \
    echo "swim ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# -----------------------------
# Switch to user
# -----------------------------
USER swim
ENV PATH="/home/swim/.local/bin:$PATH"
WORKDIR /home/swim


# -----------------------------
# Install Boost (fixed path)
# -----------------------------
RUN curl -L https://archives.boost.io/release/1.91.0/source/boost_1_91_0.tar.gz \
    -o boost.tar.gz && \
    tar -xzf boost.tar.gz && \
    rm boost.tar.gz && \
    mv boost_* boost

ENV BOOST_ROOT=/home/swim/boost

# -----------------------------
# Python dependencies
# -----------------------------
RUN python -m pip install --no-cache-dir \
    numpy scipy matplotlib scikit-learn joblib numba \
    jupyterlab cffi mpi4py cobaya cython pyside6 candl-like

# -----------------------------
# Clone SWIM
# -----------------------------
RUN git clone https://github.com/umg-kmr/SWIM.git

WORKDIR /home/swim/SWIM

# ----------------------
# Patch compile script
# -----------------------------
RUN sed -i 's|BOOST_PATH=.*|BOOST_PATH=/home/swim/boost|' compile_SWIM.sh && \
    sed -i 's|SWIM_PATH=.*|SWIM_PATH=/home/swim/SWIM|' compile_SWIM.sh

# -----------------------------
# Compile SWIM
# -----------------------------
#RUN chmod +x compile_SWIM.sh && \
#    ./compile_SWIM.sh

# -----------------------------
# Default shell
# -----------------------------
WORKDIR /home/swim/SWIM
ENV SHELL=/bin/bash
RUN echo 'source /etc/bash.bashrc' >> /home/swim/.bashrc && \
    echo 'export PS1="\\u@\\h:\\w\\$ "' >> /home/swim/.bashrc
CMD ["bash"]
