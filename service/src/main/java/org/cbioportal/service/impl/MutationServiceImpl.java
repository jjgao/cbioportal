package org.cbioportal.service.impl;

import org.cbioportal.model.*;
import org.cbioportal.model.meta.MutationMeta;
import org.cbioportal.persistence.MutationRepository;
import org.cbioportal.service.MolecularProfileService;
import org.cbioportal.service.MutationService;
import org.cbioportal.service.exception.MolecularProfileNotFoundException;
import org.cbioportal.service.util.ChromosomeCalculator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.List;

@Service
public class MutationServiceImpl implements MutationService {

    @Autowired
    private MutationRepository mutationRepository;
    @Autowired
    private MolecularProfileService molecularProfileService;
    @Autowired
    private ChromosomeCalculator chromosomeCalculator;

    @Override
    public List<Mutation> getMutationsInMolecularProfileBySampleListId(String molecularProfileId, String sampleListId,
                                                                       List<Integer> entrezGeneIds, Boolean snpOnly,
                                                                       String projection, Integer pageSize,
                                                                       Integer pageNumber, String sortBy,
                                                                       String direction)
        throws MolecularProfileNotFoundException {

        validateMolecularProfile(molecularProfileId);

        List<Mutation> mutationList = mutationRepository.getMutationsInMolecularProfileBySampleListId(molecularProfileId,
            sampleListId, entrezGeneIds, snpOnly, projection, pageSize, pageNumber, sortBy, direction);

        mutationList.forEach(mutation -> chromosomeCalculator.setChromosome(mutation.getGene()));
        return mutationList;
    }

    @Override
    public MutationMeta getMetaMutationsInMolecularProfileBySampleListId(String molecularProfileId, String sampleListId,
                                                                         List<Integer> entrezGeneIds)
        throws MolecularProfileNotFoundException {

        validateMolecularProfile(molecularProfileId);

        return mutationRepository.getMetaMutationsInMolecularProfileBySampleListId(molecularProfileId, sampleListId,
            entrezGeneIds);
    }

    @Override
    public List<Mutation> getMutationsInMultipleMolecularProfiles(List<String> molecularProfileIds, 
                                                                  List<String> sampleIds, List<Integer> entrezGeneIds, 
                                                                  String projection, Integer pageSize, 
                                                                  Integer pageNumber, String sortBy, String direction) {

        List<Mutation> mutationList = mutationRepository.getMutationsInMultipleMolecularProfiles(molecularProfileIds,
            sampleIds, entrezGeneIds, projection, pageSize, pageNumber, sortBy, direction);

        mutationList.forEach(mutation -> chromosomeCalculator.setChromosome(mutation.getGene()));
        return mutationList;
    }

    @Override
    public MutationMeta getMetaMutationsInMultipleMolecularProfiles(List<String> molecularProfileIds, 
                                                                    List<String> sampleIds, 
                                                                    List<Integer> entrezGeneIds) {

        return mutationRepository.getMetaMutationsInMultipleMolecularProfiles(molecularProfileIds, sampleIds,
            entrezGeneIds);
    }

    @Override
    public List<Mutation> fetchMutationsInMolecularProfile(String molecularProfileId, List<String> sampleIds,
                                                           List<Integer> entrezGeneIds, Boolean snpOnly,
                                                           String projection, Integer pageSize, Integer pageNumber,
                                                           String sortBy, String direction)
        throws MolecularProfileNotFoundException {

        validateMolecularProfile(molecularProfileId);

        List<Mutation> mutationList = mutationRepository.fetchMutationsInMolecularProfile(molecularProfileId, sampleIds,
            entrezGeneIds, snpOnly, projection, pageSize, pageNumber, sortBy, direction);

        mutationList.forEach(mutation -> chromosomeCalculator.setChromosome(mutation.getGene()));
        return mutationList;
    }

    @Override
    public MutationMeta fetchMetaMutationsInMolecularProfile(String molecularProfileId, List<String> sampleIds,
                                                             List<Integer> entrezGeneIds)
        throws MolecularProfileNotFoundException {

        validateMolecularProfile(molecularProfileId);

        return mutationRepository.fetchMetaMutationsInMolecularProfile(molecularProfileId, sampleIds, entrezGeneIds);
    }

    @Override
    public List<MutationCountByGene> getSampleCountByEntrezGeneIdsAndSampleIds(String molecularProfileId,
                                                                               List<String> sampleIds,
                                                                               List<Integer> entrezGeneIds)
        throws MolecularProfileNotFoundException {

        validateMolecularProfile(molecularProfileId);

        return mutationRepository.getSampleCountByEntrezGeneIdsAndSampleIds(molecularProfileId, sampleIds, 
            entrezGeneIds);
    }

    @Override
    public List<MutationCountByGene> getPatientCountByEntrezGeneIdsAndSampleIds(String molecularProfileId,
                                                                                List<String> patientIds,
                                                                                List<Integer> entrezGeneIds)
        throws MolecularProfileNotFoundException {

        validateMolecularProfile(molecularProfileId);

        return mutationRepository.getPatientCountByEntrezGeneIdsAndSampleIds(molecularProfileId, patientIds, 
            entrezGeneIds);
    }

    @Override
    public List<MutationCount> getMutationCountsInMolecularProfileBySampleListId(String molecularProfileId,
                                                                                 String sampleListId)
        throws MolecularProfileNotFoundException {

        validateMolecularProfile(molecularProfileId);

        return mutationRepository.getMutationCountsInMolecularProfileBySampleListId(molecularProfileId, sampleListId);
    }

    @Override
    public List<MutationCount> fetchMutationCountsInMolecularProfile(String molecularProfileId, List<String> sampleIds)
        throws MolecularProfileNotFoundException {

        validateMolecularProfile(molecularProfileId);

        return mutationRepository.fetchMutationCountsInMolecularProfile(molecularProfileId, sampleIds);
    }

    @Override
    public List<MutationCountByPosition> fetchMutationCountsByPosition(List<Integer> entrezGeneIds,
                                                                       List<Integer> proteinPosStarts,
                                                                       List<Integer> proteinPosEnds) {

        List<MutationCountByPosition> mutationCountByPositionList = new ArrayList<>();
        for (int i = 0; i < entrezGeneIds.size(); i++) {

            MutationCountByPosition mutationCountByPosition = mutationRepository.getMutationCountByPosition(
                entrezGeneIds.get(i), proteinPosStarts.get(i), proteinPosEnds.get(i));
            mutationCountByPositionList.add(mutationCountByPosition);
        }

        return mutationCountByPositionList;
    }

    private void validateMolecularProfile(String molecularProfileId) throws MolecularProfileNotFoundException {

        MolecularProfile molecularProfile = molecularProfileService.getMolecularProfile(molecularProfileId);

        if (!molecularProfile.getMolecularAlterationType()
            .equals(MolecularProfile.MolecularAlterationType.MUTATION_EXTENDED)) {

            throw new MolecularProfileNotFoundException(molecularProfileId);
        }
    }
}
