include { logColours           } from '../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion           } from '../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def getDBFlag(db_list, db_name, value_for_all) {
    if (db_list.contains(value_for_all)) {
        return true
    } else if (db_list.contains(db_name)) {
        return true
    } else {
        return false
    }
}


def checkDBVersion(version_file, db_version, db_name) {
    if (!version_file.exists()) {
        error("Version file for DB ${db_name} not found at: ${version_file}")
    }
    version_file.withReader { r ->
        def version = r.readLine()
        if (version != db_version) {
            error("Version for DB ${db_name} did not match expected version: ${db_version}. Version found: ${version}")
        }
    }
}