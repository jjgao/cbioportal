/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mskcc.cbio.portal.hotspots;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.mskcc.cbio.portal.dao.DaoException;
import org.mskcc.cbio.portal.dao.DaoGeneOptimized;
import org.mskcc.cbio.portal.dao.DaoPdbUniprotResidueMapping;
import org.mskcc.cbio.portal.dao.DaoProteinContactMap;
import org.mskcc.cbio.portal.dao.DaoUniProtIdMapping;
import org.mskcc.cbio.portal.model.CanonicalGene;
import org.mskcc.cbio.portal.model.PdbUniprotAlignment;
import org.mskcc.cbio.portal.model.PdbUniprotResidueMapping;

/**
 *
 * @author jgao
 */
public class CalcDistanceBetweenHotspots {
    
    
    public static void main(String[] args) throws FileNotFoundException, IOException {
        
        String input = "/Users/jgao/projects/collab/indel-hotspot-3d/esr1_regions";
        FileReader reader =  new FileReader(input);
        BufferedReader in = new BufferedReader(reader);
        
        String output = input + ".contact.txt";
        FileWriter writer = new FileWriter(output);
        BufferedWriter out = new BufferedWriter(writer);
        
        in.readLine(); // skip the first line
        for (String line = in.readLine(); line != null; line = in.readLine()) {
            StringBuilder sb = new StringBuilder(line);
            
            sb.append(line);
            
            String[] parts = line.split("\t");
            
            sb.append("\t");
            boolean contact=false;
            try {
                contact = doesContact(parts[0], Integer.parseInt(parts[1]), Integer.parseInt(parts[2]), Integer.parseInt(parts[3]), Integer.parseInt(parts[4]));
            } catch (Exception ex) {
                
            }
            sb.append(Boolean.toString(contact));
                
            sb.append("\t");
            double dist=Double.MAX_VALUE;
            try {
                dist = getDistance(parts[0], Integer.parseInt(parts[1]), Integer.parseInt(parts[2]), Integer.parseInt(parts[3]), Integer.parseInt(parts[4]));
            } catch (Exception ex) {
            }
            sb.append(Double.toString(dist));            
            
            out.append(sb);
            out.newLine();
            System.out.println(sb);
        }
        
        in.close();
        out.close();
        
    }
    
    private static Map<String, List<PdbUniprotAlignment>> alignments = new HashMap<String, List<PdbUniprotAlignment>>();
    private static List<PdbUniprotAlignment> getAlignments(String hugo) throws Exception {
        List<PdbUniprotAlignment> alns = alignments.get(hugo);
        if (alns==null) {
            alns = new ArrayList<>();
            alignments.put(hugo, alns);
            
            CanonicalGene gene = DaoGeneOptimized.getInstance().getGene(hugo);
            List<String> uniprotAccs = DaoUniProtIdMapping.mapFromEntrezGeneIdToUniprotAccession(gene.getEntrezGeneId());
            for (String uniprotAcc : uniprotAccs) {
                String uniprotId = DaoUniProtIdMapping.mapFromUniprotAccessionToUniprotId(uniprotAcc);
                alns.addAll(DaoPdbUniprotResidueMapping.getAlignments(uniprotId));
            }
        }
        return alns;
    }
    
    private static double getDistance(String hugo, int start1, int end1, int start2, int end2) throws Exception {
        double shortestDistance = Double.MAX_VALUE;
        List<PdbUniprotAlignment> alns = getAlignments(hugo);
        for (PdbUniprotAlignment aln : alns) {

            Collection<PdbUniprotResidueMapping> residueMapping1 = getResidueMapping(aln, start1, end1);
            if (residueMapping1.isEmpty()) {
                continue;
            }

            Collection<PdbUniprotResidueMapping> residueMapping2 = getResidueMapping(aln, start2, end2);
            if (residueMapping2.isEmpty()) {
                continue;
            }

            Chain chain = getChain(aln.getPdbId(), aln.getChain());

            List<Group> residues1 = new ArrayList<Group>();
            for (PdbUniprotResidueMapping m : residueMapping1) {
                ResidueNumber rn = new ResidueNumber();
                rn.setChainId(aln.getChain());
                rn.setSeqNum(m.getPdbPos());
                if (m.getPdbInsertionCode()!=null &&  m.getPdbInsertionCode()!="") {
                    rn.setInsCode(m.getPdbInsertionCode().charAt(0));
                }
                Group group = chain.getGroupByPDB(rn);
                residues1.add(group);
            }
            List<Group> residues2 = new ArrayList<Group>();
            for (PdbUniprotResidueMapping m : residueMapping2) {
                ResidueNumber rn = new ResidueNumber();
                rn.setChainId(aln.getChain());
                rn.setSeqNum(m.getPdbPos());
                if (m.getPdbInsertionCode()!=null &&  m.getPdbInsertionCode()!="") {
                    rn.setInsCode(m.getPdbInsertionCode().charAt(0));
                }
                Group group = chain.getGroupByPDB(rn);
                residues2.add(group);
            }


            for (Group r1 : residues1) {
                Atom a1 = r1.getAtom("CA");
                for (Group r2 : residues2) {
                    Atom a2 = r2.getAtom("CA");
                    double dist = Calc.getDistance(a1, a2);
                    if (shortestDistance > dist) {
                        shortestDistance = dist;
                    }
                }
            }

        }
        return shortestDistance;
    }
        
    
    private static AtomCache atomCache = getAtomCache("/Users/jgao/projects/cbio-portal-data/reference-data/pdb-cache");
    private static Map<String, Structure> strucs = new HashMap<String, Structure>();
    private static Structure getStructure(String pdbId) throws Exception {
        Structure struc = strucs.get(pdbId);
        if (struc==null) {
            struc = atomCache.getStructure(pdbId);
            strucs.put(pdbId, struc);
        }
        return struc;
    }
    
    private static Map<String, Chain> chains = new HashMap<String, Chain>();
    private static Chain getChain(String pdbId, String chainId) throws Exception {
        String key = pdbId + "_" + chainId;
        Chain chain = chains.get(key);
        if (chain==null) {
            Structure struc = getStructure(pdbId);
            chain = struc.getChainByPDB(chainId);
            chains.put(key, chain);
        }
        return chain;
    }
    
    private static AtomCache getAtomCache(String dirCache) {
        AtomCache atomCache = new AtomCache(dirCache, true);
        FileParsingParameters params = new FileParsingParameters();
//        params.setLoadChemCompInfo(true);
        params.setAlignSeqRes(true);
        params.setParseSecStruc(false);
        params.setUpdateRemediatedFiles(false);
        atomCache.setFileParsingParams(params);
        atomCache.setAutoFetch(true);
        return atomCache;
    }
    
    
    private static boolean doesContact(String hugo, int start1, int end1, int start2, int end2) throws Exception {
        List<PdbUniprotAlignment> alns = getAlignments(hugo);

        for (PdbUniprotAlignment aln : alns) {
            Collection<PdbUniprotResidueMapping> residueMapping1 = getResidueMapping(aln, start1, end1);
            if (residueMapping1.isEmpty()) {
                continue;
            }

            Collection<PdbUniprotResidueMapping> residueMapping2 = getResidueMapping(aln, start2, end2);
            if (residueMapping2.isEmpty()) {
                continue;
            }

            Set<Integer> residues1 = new HashSet<Integer>();
            for (PdbUniprotResidueMapping m : residueMapping1) {
                residues1.add(m.getPdbPos());
            }
            Set<Integer> residues2 = new HashSet<Integer>();
            for (PdbUniprotResidueMapping m : residueMapping2) {
                residues2.add(m.getPdbPos());
            }

            Set<Integer> residues = new HashSet<Integer>(residues1);
            residues.addAll(residues2);
            Map<Integer, Set<Integer>> contact = DaoProteinContactMap.getProteinContactMap(aln.getPdbId(), aln.getChain(), residues, -1, -1, -1);

            for (Integer r1 : residues1) {
               if (contact.containsKey(r1) && !Collections.disjoint(contact.get(r1), residues2)) {
                   return true;
               }
            }

            for (Integer r2 : residues2) {
               if (contact.containsKey(r2) && !Collections.disjoint(contact.get(r2), residues1)) {
                   return true;
               }
            }
        }
        return false;
    }
    
    private static Map<String,Collection<PdbUniprotResidueMapping>> residueMappings = new HashMap<String,Collection<PdbUniprotResidueMapping>>();
    private static Collection<PdbUniprotResidueMapping> getResidueMapping(PdbUniprotAlignment aln, int start, int end) throws DaoException {
        String key = "_"+aln.getAlignmentId()+"_"+start+"_"+end;
        Collection<PdbUniprotResidueMapping> rm = residueMappings.get(key);
        if (rm==null) {
            Set<Integer> uniprotPositions = new HashSet<Integer>();
            for (int i=start; i<=end; i++) {
                uniprotPositions.add(i);
            }
            rm = DaoPdbUniprotResidueMapping.mapToPdbResidues(aln.getAlignmentId(), uniprotPositions).values();
            residueMappings.put(key, rm);
        }
        
        return rm;
    }
}
