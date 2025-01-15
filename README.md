# EcologyUI

Frontend code for ecology scripts project.

## Getting Started

### Backend

To get started, cd into api/ and execute `run.sh`

### Frontend

With the backend running, open a separate terminal and cd into frontend/

To run in a browser, execute `npm run dev`

To run as a standalone app, execute `npx tauri dev`

## Tech Stack

The UI is written in React and uses Vite as a build tool. A separate backend process supervises the analyses and is written in Rust using Tauri, spawning R sessions as needed to run the analyses. 

Rmd scripts are called within each respective R session which interface with functions provided in the EcologyCore package to generate the final report. 

Reports can be generated in markdown, HTML, or Word format and can be saved as PDFs. RData files are also outputted as a non-destructive way to save results for further analysis.

Each part of the stack can be used independently, so the EcologyCore package can be used in R without any of the graphical UX if desired.

## Design Philosophy

The idea is to build a modern web/standalone application to simplify running and interpreting ecological analyses. The frontend provides a GUI to select the necessary input data and select from a list of analyses. This invokes a function in the Rust backend, using IPC, which acts as an interop layer to pass the data to an R process. The R process generates HTML output from an Rmd script, which calls the relevant functions in the EcologyCore package in R, and presents the outputted figures alongside descriptive text which provides the necessary background to interpret the figures. Moreover, an RData file is generated as a way of retaining the results of the analysis in a non-destructive way, which can then be handled by the user separately if they wish. 