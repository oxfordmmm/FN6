ARG APP_NAME=fn6

################################################################################
# Create a stage for building the application.
# Notably we need python3.12 for the build process (due to the python bindings)
################################################################################

FROM python:3.12.13-trixie AS build
ARG APP_NAME
WORKDIR /app

# Install Rust and Cargo.
RUN curl -sSL https://sh.rustup.rs | sh -s -- -y --default-toolchain stable && \
    . $HOME/.cargo/env


# Build the application.
RUN --mount=type=bind,source=src,target=src \
    --mount=type=bind,source=Cargo.toml,target=Cargo.toml \
    --mount=type=bind,source=Cargo.lock,target=Cargo.lock \
    --mount=type=bind,source=README.md,target=README.md \
    --mount=type=cache,target=/app/target/ \
    --mount=type=cache,target=/usr/local/cargo/git/db \
    --mount=type=cache,target=/usr/local/cargo/registry/ \
    . $HOME/.cargo/env && \
    cargo build --locked --release && \
    cp ./target/release/$APP_NAME /bin/fn6

# For some reason `rev` is not a part of the util-linux package in >ubuntu:questing, so pin here for now
# This is needed as part of the pipeline, but should be inconsequential to pin otherwise
FROM ubuntu:questing AS final

# Copy the executable from the "build" stage.
COPY --from=build /bin/fn6 /bin/

RUN apt update && apt install -y curl jq pigz && rm -rf /var/lib/apt/lists/*


# What the container should run when it is started.
CMD ["/bin/fn6"]