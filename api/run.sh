#!/bin/bash

# Directory to watch
WATCH_DIR="${1:-.}"

# Check if the directory exists
if [[ ! -d "$WATCH_DIR" ]]; then
    echo "Error: Directory $WATCH_DIR does not exist."
    exit 1
fi

# Store the initial state of files
declare -A file_timestamps

# PID of the currently running R session
R_SESSION_PID=""

# Function to get file modification times
update_timestamps() {
    for file in "$WATCH_DIR"/*.R; do
        if [[ -f "$file" ]]; then
            file_timestamps["$file"]=$(stat -c %Y "$file")
        fi
    done
}

# Function to check for changes
check_for_changes() {
    for file in "$WATCH_DIR"/*.R; do
        if [[ -f "$file" ]]; then
            new_time=$(stat -c %Y "$file")
            if [[ "${file_timestamps["$file"]}" != "$new_time" ]]; then
                file_timestamps["$file"]=$new_time
                return 0
            fi
        fi
    done
    return 1
}

# Function to stop the previous R session
stop_previous_r_session() {
    if [[ -n "$R_SESSION_PID" ]] && kill -0 "$R_SESSION_PID" 2>/dev/null; then
        echo "Stopping previous R session (PID: $R_SESSION_PID)..."
        kill "$R_SESSION_PID"
        wait "$R_SESSION_PID" 2>/dev/null
        echo "Previous R session stopped."
    fi
}

# Function to execute R scripts
execute_r_scripts() {
    echo "Executing R scripts in $WATCH_DIR..."
    Rscript api.R &
    R_SESSION_PID=$!
    # echo "Executing R scripts in $WATCH_DIR..."
    # for script in "$WATCH_DIR"/*.R; do
    #     if [[ -f "$script" ]]; then
    #         echo "Running $script..."
    #         Rscript "$script"
    #     fi
    # done
}

# Initialize timestamps
update_timestamps

# Initial execution
execute_r_scripts

echo "Watching directory: $WATCH_DIR for changes to R scripts..."

# Watch for changes in a loop
while true; do
    if check_for_changes; then
        echo "Change detected in R scripts."
        execute_r_scripts
    fi
    sleep 1  # Check every 1 second
done