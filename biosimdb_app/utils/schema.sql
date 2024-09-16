/* Create the Tables */
/* These need to be created before you can store and retrieve data. 
Create a file with the SQL commands needed to create empty tables: */

DROP TABLE IF EXISTS creator;
DROP TABLE IF EXISTS project;
DROP TABLE IF EXISTS simulation;
DROP TABLE IF EXISTS topology;
DROP TABLE IF EXISTS trajectory;
DROP TABLE IF EXISTS aiida_archive;

CREATE TABLE creator (
  creator_ID INTEGER PRIMARY KEY AUTOINCREMENT,
  creator VARCHAR(140) NOT NULL,
  email VARCHAR(140) UNIQUE NOT NULL
);

create TABLE project (
  project_ID INTEGER PRIMARY KEY AUTOINCREMENT,
  creator_ID INTEGER NOT NULL,
  `datetime added` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  authors VARCHAR(500) NOT NULL, /* string of author list */
  title VARCHAR(140) NOT NULL,
  abstract VARCHAR(2000) NOT NULL,
  citations VARCHAR(500), /* string of DOI list */
  FOREIGN KEY (creator_ID) REFERENCES creator (creator_ID)
);

CREATE TABLE simulation (
  simulation_ID INTEGER PRIMARY KEY AUTOINCREMENT,
  project_ID INTEGER NOT NULL,
  temperature FLOAT,
  pressure FLOAT,
  `number of frames` FLOAT,
  timestep FLOAT,
  thermostat VARCHAR(140),
  barostat VARCHAR(140),
  software VARCHAR(140),
  forcefields VARCHAR(500),
  `experimental structures` VARCHAR(500),
  `simulation directory` VARCHAR(140),
  FOREIGN KEY (project_ID) REFERENCES project (project_ID)
);

CREATE TABLE topology (
  topology_ID INTEGER PRIMARY KEY AUTOINCREMENT,
  simulation_ID INTEGER NOT NULL,
  `atoms list` TEXT, /* dict */
  `masses list` TEXT, /* dict */
  `molecules list` TEXT, /* dict */
  `charges list` TEXT, /* dict */
  `system charge` FLOAT,
  FOREIGN KEY (simulation_ID) REFERENCES simulation (simulation_ID)
);

CREATE TABLE trajectory (
  trajectory_ID INTEGER PRIMARY KEY AUTOINCREMENT,
  simulation_ID INTEGER NOT NULL,
  `number of atoms`s INT,
  `number of frames` INT,
  length FLOAT,
  `time between snapshots` FLOAT,
  `single frame coordinates` TEXT, /* dict */
  `single frame dimensions` TEXT, /* dict */
  FOREIGN KEY (simulation_ID) REFERENCES simulation (simulation_ID)
);

CREATE TABLE aiida_archive (
  aiida_archive_ID INTEGER PRIMARY KEY AUTOINCREMENT,
  simulation_ID INTEGER NOT NULL,
  `aiida version` VARCHAR(140),
  `aiida plugin` VARCHAR(140),
  FOREIGN KEY (simulation_ID) REFERENCES simulation (simulation_ID)
);