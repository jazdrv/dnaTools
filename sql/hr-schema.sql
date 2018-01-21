CREATE TABLE IF NOT EXISTS `FeatureType` (
  `id` smallint(6) NOT NULL AUTO_INCREMENT,
  `description` varchar(64) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `id` (`id`,`description`)
);

CREATE TABLE IF NOT EXISTS `Contig` (
  `id` smallint(6) NOT NULL AUTO_INCREMENT,
  `build` char(4) NOT NULL,
  `description` varchar(64) NOT NULL,
  `length` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `build_desc_idx` (`build`,`description`)
);

CREATE TABLE IF NOT EXISTS `Feature` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `description` varchar(64) NOT NULL,
  `ref` varchar(255) NOT NULL,
  `alt` varchar(255) NOT NULL,
  `refDerived` tinyint(1) NOT NULL,
  `public` tinyint(1) NOT NULL,
  `combBED` tinyint(1) NOT NULL,
  `featureTypeId` smallint(6) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `desc` (`description`),
  KEY `combBED` (`combBED`),
  KEY `featureTypeId` (`featureTypeId`)
);

CREATE TABLE IF NOT EXISTS `FeatureAlias` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `featureId` int(11) NOT NULL,
  `name` varchar(64) NOT NULL,
  `reference` varchar(64) NOT NULL,
  `year` smallint(6) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`),
  KEY `featureId` (`featureId`)
);

CREATE TABLE IF NOT EXISTS `FeatureLocation` (
  `id` bigint(20) NOT NULL AUTO_INCREMENT,
  `featureID` int(11) NOT NULL,
  `contigID` smallint(6) NOT NULL,
  `start` int(11) NOT NULL,
  `stop` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `featureID_2` (`featureID`,`contigID`),
  UNIQUE KEY `featureID` (`featureID`),
  KEY `contigID` (`contigID`)
);

CREATE TABLE IF NOT EXISTS `TreeLeaf` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `parentId` int(11) DEFAULT NULL,
  `description` varchar(64) NOT NULL,
  `tmrca` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `description` (`description`),
  KEY `parentId` (`parentId`)
);

CREATE TABLE IF NOT EXISTS `TreeEdge` (
  `ancestorId` int(11) NOT NULL,
  `descendantId` int(11) NOT NULL,
  UNIQUE KEY `primary_idx` (`ancestorId`,`descendantId`),
  KEY `ancestorId` (`ancestorId`),
  KEY `descendantId` (`descendantId`)
);

CREATE TABLE IF NOT EXISTS `TreeLeafFeatureEdge` (
  `leafId` int(11) NOT NULL,
  `featureId` int(11) NOT NULL,
  UNIQUE KEY `primary_idx` (`leafId`,`featureId`),
  KEY `leafId` (`leafId`),
  KEY `featureId` (`featureId`)
);

CREATE TABLE IF NOT EXISTS `FeatureOccurrence` (
  `featureId` int(11) NOT NULL,
  `treeLeafId` int(11) NOT NULL,
  `derivedCnt` int(11) NOT NULL,
  `ancestralCnt` int(11) NOT NULL,
  UNIQUE KEY `featureId` (`featureId`,`treeLeafId`)
);

CREATE TABLE IF NOT EXISTS `Person` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `firstName` varchar(64) CHARACTER SET utf16 DEFAULT NULL,
  `middleName` varchar(64) CHARACTER SET utf16 DEFAULT NULL,
  `surname` varchar(64) CHARACTER SET utf16 DEFAULT NULL,
  `maidenName` varchar(64) CHARACTER SET utf16 DEFAULT NULL,
  `yHaplogroupId` int(11) DEFAULT NULL,
  `mtHaplogroupId` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `mtHaplogroupId` (`mtHaplogroupId`),
  KEY `yHaplogroupId` (`yHaplogroupId`)
);

CREATE TABLE IF NOT EXISTS `Lab` (
  `id` smallint(6) NOT NULL AUTO_INCREMENT,
  `description` varchar(64) NOT NULL,
  PRIMARY KEY (`id`)
);

CREATE TABLE IF NOT EXISTS `TestType` (
  `id` smallint(6) NOT NULL AUTO_INCREMENT,
  `description` varchar(64) NOT NULL,
  `isNGS` tinyint(1) NOT NULL,
  `tagNm` varchar(32) NOT NULL,
  `priority` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `description_idx` (`description`)
);

CREATE TABLE IF NOT EXISTS `LabTestEdge` (
  `labId` smallint(6) NOT NULL,
  `testTypeId` smallint(6) NOT NULL,
  UNIQUE KEY `lab_test_indx` (`labId`,`testTypeId`),
  KEY `labId` (`labId`),
  KEY `testTypeId` (`testTypeId`)
);

CREATE TABLE IF NOT EXISTS `Kit` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `labId` smallint(6) NOT NULL,
  `description` varchar(64) NOT NULL,
  `testId` smallint(6) NOT NULL,
  `callable` int(11) NOT NULL,
  `no_coverage` int(11) NOT NULL,
  `low_coverage` int(11) NOT NULL,
  `poor_coverage` int(11) NOT NULL,
  `combBED` int(11) NOT NULL,
  `pRegion` int(11) NOT NULL,
  `personId` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `lab_id_index` (`labId`,`description`),
  KEY `personId` (`personId`),
  KEY `labId` (`labId`),
  KEY `testId` (`testId`)
);

CREATE TABLE IF NOT EXISTS `KitFeature` (
  `kitId` int(11) NOT NULL,
  `featureId` int(11) NOT NULL,
  `allele` varchar(255) NOT NULL,
  `depth` tinyint(4) NOT NULL,
  `gq` tinyint(4) NOT NULL,
  UNIQUE KEY `kitfeature_idx` (`featureId`,`kitId`,`allele`),
  KEY `kitId` (`kitId`)
);

CREATE TABLE IF NOT EXISTS `UploadLog` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `contact` varchar(256) NOT NULL,
  `kitID` varchar(64) NOT NULL,
  `surname` varchar(64) CHARACTER SET utf8 DEFAULT NULL,
  `country` varchar(256) CHARACTER SET utf8 DEFAULT NULL,
  `birthYr` int(11) DEFAULT NULL,
  `labID` smallint(6) NOT NULL,
  `testTypeID` smallint(6) NOT NULL,
  `buildNm` char(3) NOT NULL,
  `importDt` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `fileNm` varchar(256) NOT NULL,
  `origFileNm` varchar(256) DEFAULT NULL,
  `otherInfo` varchar(2048) CHARACTER SET utf8 DEFAULT NULL,
  `normalOrig` varchar(256) CHARACTER SET utf8 DEFAULT NULL,
  `lat` decimal(11,8) DEFAULT NULL,
  `lng` decimal(11,8) DEFAULT NULL,
  `policyVer` varchar(64) DEFAULT NULL,
  `accessToken` varchar(1024) DEFAULT NULL,
  `updated` timestamp NULL DEFAULT NULL ON UPDATE CURRENT_TIMESTAMP,
  `hold` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `kit_lab_test_build_idx` (`kitID`,`labID`,`testTypeID`,`buildNm`),
  KEY `contact` (`contact`),
  KEY `labInd` (`labID`),
  KEY `testTypeInd` (`testTypeID`)
);

