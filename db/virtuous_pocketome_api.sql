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
-- Database: `virtuous_pocketome_api`
--

-- --------------------------------------------------------

--
-- Table structure for table `pocketome_jobs`
--

DROP TABLE IF EXISTS `pocketome_jobs`;
CREATE TABLE IF NOT EXISTS `pocketome_jobs` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `submitted` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `started` timestamp NULL DEFAULT NULL,
  `completed` timestamp NULL DEFAULT NULL,
  `user` varchar(255) CHARACTER SET utf8mb4 NOT NULL,
  `complex_pdb` varchar(255) CHARACTER SET utf8mb4 NOT NULL,
  `chain_protein` varchar(1) CHARACTER SET utf8mb4 NOT NULL,
  `chain_ligand` varchar(1) CHARACTER SET utf8mb4 NOT NULL,
  `complex_xtc` varchar(255) CHARACTER SET utf8mb4 DEFAULT NULL,
  `k_flag` int(11) DEFAULT NULL,
  `dist` float DEFAULT NULL,
  `sasa_threshold` float DEFAULT NULL,
  `dock_threshold` float DEFAULT NULL,
  `extensive` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
COMMIT;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
