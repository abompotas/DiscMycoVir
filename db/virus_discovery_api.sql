-- phpMyAdmin SQL Dump
-- version 5.0.2
-- https://www.phpmyadmin.net/
--
-- Host: 127.0.0.1:3306
-- Generation Time: Sep 14, 2023 at 01:53 PM
-- Server version: 5.7.31
-- PHP Version: 7.4.9

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
START TRANSACTION;
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;

--
-- Database: `virus_discovery_api`
--

-- --------------------------------------------------------

--
-- Table structure for table `virus_discovery_jobs`
--

DROP TABLE IF EXISTS `virus_discovery_jobs`;
CREATE TABLE IF NOT EXISTS `virus_discovery_jobs` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `submitted` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `started` timestamp NULL DEFAULT NULL,
  `completed` timestamp NULL DEFAULT NULL,
  `user` varchar(255) CHARACTER SET utf8mb4 NOT NULL,
  `genome` varchar(255) CHARACTER SET utf8mb4 NOT NULL,
  `paired` boolean NOT NULL DEFAULT 0,
  `adapter` varchar(255) CHARACTER SET utf8mb4 DEFAULT "NexteraPE-PE.fa:2:30:10",
  `sample_name` varchar(255) CHARACTER SET utf8mb4 DEFAULT "Sample X",
  `min_len` int(6) DEFAULT 50,
  `window` varchar(20) CHARACTER SET utf8mb4 DEFAULT "(5:20)",
  `forward_file` varchar(255) CHARACTER SET utf8mb4 DEFAULT NULL,
  `reverse_file` varchar(255) CHARACTER SET utf8mb4 DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
COMMIT;

-- Inserting a job for single-end sequencing
INSERT INTO `virus_discovery_jobs` (`submitted`, `user`, `genome`, `paired`, `adapter`, `sample_name`, `min_len`, `window`)
VALUES (CURRENT_TIMESTAMP, 'user1', 'AGTCGATCGATCGATC', 0, 'NexteraPE-PE.fa:2:30:10"', 'sample1', 50, '1:15');

-- Inserting a job for paired-end sequencing
INSERT INTO `virus_discovery_jobs` (`submitted`, `user`, `genome`, `paired`, `adapter`, `sample_name`, `min_len`, `window`)
VALUES (CURRENT_TIMESTAMP, 'user2', 'AGTCGATCGATCGATC', 1, 'NexteraPE-PE.fa:2:30:10"', 'sample2', 50, '2:16');

-- Inserting another job for single-end sequencing
INSERT INTO `virus_discovery_jobs` (`submitted`, `user`, `genome`, `paired`, `adapter`, `sample_name`, `min_len`, `window`)
VALUES (CURRENT_TIMESTAMP, 'user3', 'AGTCGATCGATCGATC', 0, 'NexteraPE-PE.fa:2:30:10"', 'sample3', 50, '3:17');

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
