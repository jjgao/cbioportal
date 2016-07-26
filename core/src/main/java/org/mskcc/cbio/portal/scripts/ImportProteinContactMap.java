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
package org.mskcc.cbio.portal.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.mskcc.cbio.portal.dao.DaoException;
import org.mskcc.cbio.portal.dao.DaoProteinContactMap;
import org.mskcc.cbio.portal.dao.MySQLbulkLoader;
import org.mskcc.cbio.portal.util.ConsoleUtil;
import org.mskcc.cbio.portal.util.FileUtil;
import org.mskcc.cbio.portal.util.ProgressMonitor;

/**
 *
 * @author jgao
 */
public class ImportProteinContactMap {
    
    public static void importData(File file, double closestAtomCutoff, double calphaCutoff) throws IOException, DaoException {
        Pattern patternRes = Pattern.compile("[A-Z]{3}:(-?[0-9]+)");
        MySQLbulkLoader.bulkLoadOn();
        FileReader reader = new FileReader(file);
        BufferedReader buf = new BufferedReader(reader);
        String line;
        while ((line = buf.readLine()) != null) {
            ProgressMonitor.incrementCurValue();
            ConsoleUtil.showProgress();

            if (!line.startsWith("#")) {
                String parts[] = line.split("\t");
                
                double distanceClosestAtoms = Double.parseDouble(parts[6]);
                if (distanceClosestAtoms>closestAtomCutoff) {
                    continue;
                }
                
                double distanceCAlpha = Double.parseDouble(parts[7]);
                if (distanceCAlpha>calphaCutoff) {
                    continue;
                }
                
                String pdbId = parts[0];
                String chainId = parts[1];
                
                // residue 1
                Matcher m = patternRes.matcher(parts[2]);
                if (!m.matches()) {
                    continue;
                }
                int res1 = Integer.parseInt(m.group(1));
                
                String atom1 = parts[3];
                
                // residue 2
                m = patternRes.matcher(parts[4]);
                if (!m.matches()) {
                    continue;
                }
                int res2 = Integer.parseInt(m.group(1));
                
                String atom2 = parts[5];
                
                double error = Double.parseDouble(parts[8]);
                
                DaoProteinContactMap.addProteinContactMap(pdbId, chainId, res1, atom1, res2, atom2, distanceClosestAtoms, distanceCAlpha, error);
            }
        }
        if (MySQLbulkLoader.isBulkLoad()) {
           MySQLbulkLoader.flushAll();
        }        
    }
    
    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.out.println("command line usage:  importPdbContactMap.pl <pdb-contact-map.txt> <closest-atom-distance-cutoff> <c-alpha-distance-cutoff>");
            System.exit(1);
        }
        DaoProteinContactMap.deleteAllRecords();
        ProgressMonitor.setConsoleMode(true);

        File file = new File(args[0]);
        double closestAtomCutoff = Double.parseDouble(args[1]);
        double calphaCutoff = Double.parseDouble(args[2]);
        System.out.println("Reading data from:  " + file.getAbsolutePath());
        int numLines = FileUtil.getNumLines(file);
        System.out.println(" --> total number of lines:  " + numLines);
        ProgressMonitor.setMaxValue(numLines);
        importData(file, closestAtomCutoff, calphaCutoff);
        ConsoleUtil.showWarnings();
        System.err.println("Done.");
    }
}
