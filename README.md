Trials Replication through Observational study across OHDSI (TROY)
=============================================================================

- Analytics use case(s): **Population-Level Estimation**
- Study type: **Clinical Application**
- Tags: **-**
- Study lead: **Kyung Won Kim**, **Seng Chan You**
- Study lead forums tag: **[Seng Chan You](https://dr-you-group.github.io/profiles/)**
- Study start date: **1 January 2022**
- Study end date: **-**
- Protocol: **[Protocol](https://github.com/ohdsi-studies/Troy/blob/master/documents/TROY_Protocol_v0.4.docx)**
- Publications: **-**
- Results explorer: **-**

We initiated ‘Trials Replication through Observational study across OHDSI (TROY)’ project to generate large population-level evidence for 15 pivotal RCTs in the real world: Type 2 diabetes mellitus, atrial arrhythmia, acute coronary syndrome, and rheumatoid arthritis. We select 15 pivotal trials and, to assess the generalizability to clinical practice, replicate the populations included in each trial as closely as possible. We evaluate external validity using a systematic three-dimensional approach: 1) proportion of patients meeting eligibility criteria among patients meeting FDA-approved indications; 2) a comparison of clinical characteristic differences between clinical trial patients and replicated patients; and 3) a comparative effectiveness analysis between target and comparator drugs.

How to run
============
This **[file](https://drive.google.com/file/d/1VPr-7BGyj9sU7qHMb488zKGfCyHarrlB/view?usp=sharing)** is a docker image that revised errors in StrategusInstantiatedModules.

When you open this file, you can see StrategusInstantiatedModules that are already installed and Strategus project will skip downloading modules step.

The usage instructions for this file are as follows:

 

1. Install Docker on your computer.

 

2. Open a Command Prompt (CMD) window and build the Docker file into an image.
   Modify the code below and enter it in CMD:
   ```
   docker load -i ./docker/file/path/file_name.tar
   ```
   Example: `docker load -i c:/Users/jc/troyStrategus.tar`

 

3. Verify the built image by entering the following code in CMD:
   ```
   docker images
   ```
   You can see the image named "troy:1.0".

 

4. Build the container by modifying the code below and entering it in CMD:
   ```
   docker run -it -p 8787:8787 --name container_name -e PASSWORD=pass_word image_name
   ```
   Example: `docker run -it -p 8787:8787 --name r_sos_fqs -e PASSWORD=abcd troy:1.0`

 

5. Open a web browser (e.g., Chrome) and navigate to `localhost:8787` to access RStudio.

 

   If prompted to enter an ID and password, use the following credentials:
   - ID: rstudio
   - PW: abcd or the password you specified after `-e PASSWORD=`

 

6. Once RStudio is launched, install your research project in the `/home/rstudio/git` directory. 

 

For example, you can clone a project repository using Terminal in Rstudio:
   ```
   git clone https://github.com/ohdsi-studies/troy.git
   ```
* The following instructions provide an example using "troy"


 

7. Open the `/home/git/troy/KeyringSetup.R` file and execute the following code to open the `.Renviron` file:
   ```
   usethis::edit_r_environ()
   ```

 

   In the `.Renviron` file, enter your GitHub Personal Access Token (GITHUB_PAT). Save the file and restart the R session. 
   
   If you don't have a token, execute the following code and enter the generated token:
   ```
   install.packages("usethis")
   library(usethis)
   create_github_token(scopes = c("(no scope)"), description = "R:GITHUB_PAT", host = "https://github.com")
   ```

 

8. In the `KeyringSetup.R` file, modify the code below to match your database (DB) settings and source the file:
   ```
   # Provide your environment-specific values ------
   dbms <- "redshift"
   connectionString <- "jdbc:redshift://your.server.goes.here:5439/your_cdm_database"
   username <- "username-goes-here"
   password = "password-goes-here"
   ```

 

9. Open the `/home/rstudio/git/troy/StrategusCodeToRun.R` file and modify the code below according to your DB:
    ```
    ##=========== START OF INPUTS ==========
    connectionDetailsReference <- "YUHS"
    workDatabaseSchema <- 'jaehyeongcho'
    cdmDatabaseSchema <- 'YUHS_CDM'
    outputLocation <- '~/output' #DO NOT TOUCH
    minCellCount <- 5
    cohortTableName <- "troy"
    ```
 

 
 

10. Execute the rest of the `StrategusCodeToRun.R` file.
11. If all code ran without errors, open the `/home/rstudio/git/troy/StrategusResultsUpload.R` file, run all lines.
12. The two files created in the following path are the analysis results: `/home/rstudio/output/results.sqlite` and `/home/rstudio/output/troyResults.zip`
 

Requirements
============

- A database in [Common Data Model version 5](https://github.com/OHDSI/CommonDataModel) in one of these platforms: SQL Server, Oracle, PostgreSQL, IBM Netezza, Apache Impala, Amazon RedShift, Google BigQuery, or Microsoft APS.
- 100 GB of free disk space

License
=======
The `TROY` package is licensed under Apache License 2.0

Development
===========
`TROY` was developed in ATLAS and R Studio.

### Development status

Under development.
