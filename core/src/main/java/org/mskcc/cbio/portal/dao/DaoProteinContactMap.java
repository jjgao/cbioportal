/** Copyright (c) 2012 Memorial Sloan-Kettering Cancer Center.
**
** This library is free software; you can redistribute it and/or modify it
** under the terms of the GNU Lesser General Public License as published
** by the Free Software Foundation; either version 2.1 of the License, or
** any later version.
**
** This library is distributed in the hope that it will be useful, but
** WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
** MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
** documentation provided hereunder is on an "as is" basis, and
** Memorial Sloan-Kettering Cancer Center 
** has no obligations to provide maintenance, support,
** updates, enhancements or modifications.  In no event shall
** Memorial Sloan-Kettering Cancer Center
** be liable to any party for direct, indirect, special,
** incidental or consequential damages, including lost profits, arising
** out of the use of this software and its documentation, even if
** Memorial Sloan-Kettering Cancer Center 
** has been advised of the possibility of such damage.  See
** the GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with this library; if not, write to the Free Software Foundation,
** Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
**/
package org.mskcc.cbio.portal.dao;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author jgao
 */
public class DaoProteinContactMap {
    private DaoProteinContactMap() {}
    
    public static int addProteinContactMap(String pdbId,
            String chain, int residue1, String atom1, int residue2, String atom2,
            double distanceClosestAtoms, double distanceCAlpha, double distanceError) {
        if (!MySQLbulkLoader.isBulkLoad()) {
            throw new IllegalStateException("only bulk load mode is supported ");
        } 
        //  write to the temp file maintained by the MySQLbulkLoader
        MySQLbulkLoader.getMySQLbulkLoader("protein_contact_map").insertRecord(
                pdbId,
                chain,
                Integer.toString(residue1),
                atom1,
                Integer.toString(residue2),
                atom2,
                Double.toString(distanceClosestAtoms),
                Double.toString(distanceCAlpha),
                Double.toString(distanceError));

        // return 1 because normal insert will return 1 if no error occurs
        return 1;
    }
    
    /**
     * Retrieve contact map for a set of residues in a PDB chain
     * @param pdbId
     * @param chain
     * @param residues
     * @param distanceErrorThreshold 
     * @return Map<residue, set <contacting residues>>
     */
    public static Map<Integer, Set<Integer>> getProteinContactMap(String pdbId, String chain,
            Collection<Integer> residues, double distanceClosestAtomsThreshold, double distanceCAlphaThreshold, double distanceErrorThreshold) throws DaoException {
        Connection con = null;
        PreparedStatement pstmt = null;
        ResultSet rs = null;
        try {
            con = JdbcUtil.getDbConnection(DaoProteinContactMap.class);
            
            String sql = "SELECT  `RESIDUE1`, `RESIDUE2` "
                    + "FROM  `protein_contact_map` "
                    + "WHERE `PDB_ID`='" + pdbId + "' "
                    + "AND `CHAIN`='" + chain + "' ";
            if (residues!=null) {
                String strResidues = StringUtils.join(residues, ",");
                sql += "AND `RESIDUE1` IN (" + strResidues + ") "
                    + "AND `RESIDUE2` IN (" + strResidues + ") ";
            }
            if (distanceClosestAtomsThreshold>0) {
                sql += "AND `DISTANCE_CLOSEST_ATOMS`<"+distanceClosestAtomsThreshold+" ";
            }
            if (distanceCAlphaThreshold>0) {
                sql += "AND `DISTANCE_C_ALPHA`<"+distanceCAlphaThreshold+" ";
            }
            if (distanceErrorThreshold>0) {
                sql += "AND `DISTANCE_ERROR`<"+distanceErrorThreshold+" ";
            }
            pstmt = con.prepareStatement(sql);
            rs = pstmt.executeQuery();
            
            Map<Integer, Set<Integer>> map = new HashMap<Integer, Set<Integer>>();
            
            while (rs.next()) {
                int res1 = rs.getInt(1);
                int res2 = rs.getInt(2);
                Set<Integer> set1 = map.get(res1);
                if (set1==null) {
                    set1 = new HashSet<Integer>();
                    map.put(res1, set1);
                }
                set1.add(res2);
                Set<Integer> set2 = map.get(res2);
                if (set2==null) {
                    set2 = new HashSet<Integer>();
                    map.put(res2, set2);
                }
                set2.add(res1);
            }
            
            return map;
        } catch (SQLException e) {
            throw new DaoException(e);
        } finally {
            JdbcUtil.closeAll(DaoProteinContactMap.class, con, pstmt, rs);
        }
    }
    
    public static void deleteAllRecords() throws DaoException {
        Connection con = null;
        PreparedStatement pstmt = null;
        ResultSet rs = null;
        try {
            con = JdbcUtil.getDbConnection(DaoProteinContactMap.class);
            pstmt = con.prepareStatement("TRUNCATE TABLE protein_contact_map");
            pstmt.executeUpdate();
        } catch (SQLException e) {
            throw new DaoException(e);
        } finally {
            JdbcUtil.closeAll(DaoProteinContactMap.class, con, pstmt, rs);
        }
    }
}
