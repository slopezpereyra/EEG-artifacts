# Use the official Debian base image
FROM debian:stable

# Update the package lists
RUN apt-get update

# Install system dependencies for your R package
RUN apt-get install -y \
    libgdal-dev gdal-bin libproj-dev proj-data proj-bin libgeos-dev \
    libcurl4-openssl-dev libssl-dev libfontconfig1-dev libxml2-dev  \ 
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev


# Clean up after installation to reduce the image size
RUN apt-get clean

# Use the official R image from Docker Hub
FROM r-base:latest

# Copy system dependencies from the builder stage
COPY --from=builder /usr/local /usr/local

RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_github('slopezpereyra/EEG-toolkit')"


# Set the working directory to your package directory
WORKDIR /package


