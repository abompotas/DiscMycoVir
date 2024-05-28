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
  `stage` tinyint(1) NOT NULL DEFAULT 0,
  `started_analysis` timestamp NULL DEFAULT NULL,
  `completed_analysis` timestamp NULL DEFAULT NULL,
  `started_trimming` timestamp NULL DEFAULT NULL,
  `completed_trimming` timestamp NULL DEFAULT NULL,
  `started_discovery` timestamp NULL DEFAULT NULL,
  `completed_discovery` timestamp NULL DEFAULT NULL,
  `user` varchar(255) CHARACTER SET utf8mb4 NOT NULL,
  `sample_name` varchar(255) CHARACTER SET utf8mb4 NOT NULL,
  `genome` varchar(255) CHARACTER SET utf8mb4 NOT NULL,
  `paired` boolean NOT NULL DEFAULT 0,
  `forward_file` varchar(255) CHARACTER SET utf8mb4 NOT NULL,
  `reverse_file` varchar(255) CHARACTER SET utf8mb4 DEFAULT NULL,
  `adapter` varchar(255) CHARACTER SET utf8mb4 DEFAULT NULL,
  `min_len` int(6) DEFAULT NULL,
  `window` varchar(20) CHARACTER SET utf8mb4 DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
COMMIT;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
