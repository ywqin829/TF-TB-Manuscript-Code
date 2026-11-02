"""Deterministic editorial restructuring; no statistical results are recomputed."""
import re
from copy import deepcopy

MOVED = {}
MAP = {3: 4, 4: 5, 5: 6, 6: 3, 7: 'S5'}

def remap(text):
    return re.sub(r'\b(Fig\. |Figure )([3-7])(?=[A-F\s.,;)])', lambda m: m[1] + str(MAP[int(m[2])]), text)

def revise_main(document):
    p = list(document.paragraphs)
    assert p[46].text.startswith('3.1 Donor-blocked')
    assert p[82].text.startswith('3.8 In vivo')
    assert p[88].text == '4. Discussion'
    MOVED['bal'] = [p[i].text for i in (83,84,85)]
    MOVED['bal'] = [
        text.replace(
            'Raw-count threshold sensitivity is reported in the Results.',
            'Raw-count threshold sensitivity is reported in the Supplementary Information.'
        )
        for text in MOVED['bal']
    ]
    MOVED['caption'] = remap(p[87].text)
    MOVED['network'] = p[99].text
    # Reuse original paragraphs and embedded images within the same DOCX package.
    def para(text, style=None):
        q = document.add_paragraph(text, style=style)
        q._p.getparent().remove(q._p)
        return q._p
    def old(i): return deepcopy(p[i]._p)
    split = p[78].text.index('In GSE31348, LM22/CIBERSORT')
    time = p[78].text[:split].strip() + ' This association supports coordinated change after time adjustment, but does not establish causation.'
    fractions = p[78].text[split:].replace(' The time-adjusted association supports coordinated change but cannot establish causation or independent biomarker value.', '')
    sections = [
        para('3.1 Reproducible transcriptional recovery during tuberculosis treatment', 'Heading 2'),
        old(47), old(48), para(p[57].text.split(' TBX21 and NFATC2')[0]), old(49), old(50),
        para('3.2 Immune-cell composition is associated with transcriptional recovery', 'Heading 2'),
        para(fractions), old(76),
        para('Both deconvolution methods indicated a modest rise in relative B-cell fraction alongside a larger fall in neutrophil fraction. Their concordance supports a composition-associated interpretation of the recovery program, while neither method independently measures cell abundance. The consequences for TF-regulon inference and the RNA score are summarized in Fig. 6 and examined below.'),
        old(80), old(81),
        para('3.3 Composition adjustment attenuates inferred TF-regulon recovery', 'Heading 2'),
        old(52), old(53), old(54), old(55), old(60),
        para('A second, cross-sectional BAL dataset (GSE326212) also localized PAX5 and MEF2B to B-lineage cells and TBX21 and RUNX3 to T/NK compartments. These maps provide cell-type context, not longitudinal evidence of intracellular recovery. Detailed BAL comparisons, including STAT4, are reported in Supplementary Fig. S5; the orthogonal network prior is reported in Supplementary Fig. S1.'),
        old(61), old(62), old(74), old(75), old(64), old(65), old(66), old(67),
        para('3.4 Limited evidence for incremental outcome information from TBRI beyond estimated B-cell composition', 'Heading 2'),
        para('The four-gene TB Resolution Index (TBRI) was selected from combinations ranked in GSE31348 as a compact, interpretable treatment-associated score, not as a uniquely optimal model. Consequently, that cohort does not provide independent score-performance estimates. GSE89403 did not select the TBRI genes or signs and was used to examine outcome association; it had also contributed to TF screening.'),
        old(69), old(70), old(77), para(time),
        para('The absence of consistent early improvement over estimated B-cell fraction limits evidence for incremental outcome information. It does not establish statistical equivalence or exclude a clinically useful effect in a larger prospective cohort.'),
        old(71), old(72),
    ]
    for element in sections: p[45]._p.getparent().insert(p[88]._p.getparent().index(p[88]._p), element)
    for q in p[46:88]: q._p.getparent().remove(q._p)
    p[3].text = ('Objective: To examine how immune-cell composition affects the interpretation of blood transcriptional recovery during tuberculosis treatment. Methods: Three longitudinal cohorts were analysed using participant-aware models, fixed gene-set testing, two deconvolution methods, and composition-adjusted TF-regulon and RNA-score analyses. Single-cell datasets provided cross-sectional localization. Results: In GSE31348, 162/171 discovery-defined recovery genes increased and 156/163 measurable disease-associated genes decreased during treatment. CIBERSORT and quanTIseq both estimated rising relative B-cell fractions and falling neutrophil fractions. Composition adjustment attenuated PAX5, TBX21 and MAF treatment effects; no screened TF retained all cross-cohort criteria after adjustment. The exploratory four-gene TBRI tracked estimated B-cell fraction even after categorical treatment-time adjustment (within-participant coefficient 1.816, 95% CI 1.548-2.084). It showed limited evidence of incremental early outcome information beyond estimated B-cell fraction. Conclusions: Treatment-associated blood transcriptional recovery is reproducible, but its interpretation requires accounting for immune-cell composition when inferring TF activity and evaluating RNA biomarkers. Estimates from shared transcriptomes do not establish absolute B-cell reconstitution or intracellular regulatory mechanisms.')
    p[4].text = 'Keywords: tuberculosis; transcriptional recovery; immune-cell composition; transcription factors; B cells; RNA biomarkers'
    p[13].text = ('We therefore asked how accounting for immune-cell composition changes the interpretation of treatment-associated blood transcriptional recovery. Using participants as the unit of inference, we first tested discovery-defined recovery genes in a separate longitudinal cohort, then assessed estimated cell-fraction changes and their association with the recovery program. We examined the consequences for TF-regulon inference and used the exploratory TBRI as an example of RNA biomarker interpretation. Single-cell localization supplied complementary cell-type context. These analyses distinguish reproducible expression patterns from evidence for intracellular regulation or clinical prediction.')
    # Preserve the original contextual paragraphs for the supplementary document, then
    # replace the main-text Discussion as a coherent interpretive section.  The prior
    # implementation edited fixed paragraph indices here, which silently overwrote the
    # Discussion generated by 08_build_manuscript.py with result-like fragments.
    MOVED['context'] = [p[i].text for i in range(93,98)]
    for i in range(93,100): p[i]._p.getparent().remove(p[i]._p)

    final_discussion = [
        'Earlier longitudinal studies established that blood transcriptional profiles change rapidly during effective tuberculosis therapy, with an early decline in inflammatory and interferon-associated signals followed by a slower return of lymphocyte-related programs [9,53]. The question raised by those observations is no longer simply whether the blood transcriptome responds to treatment, but what biological process a returning signature represents. Our reanalysis suggests that apparent transcriptional recovery in whole blood should be understood as an integrated change in the sampled immune system, rather than automatically as restoration of gene regulation within individual cells. This distinction matters because a reproducible change in expression can be informative about disease state while remaining ambiguous about the cellular process that generated it.',
            'The opposing trajectories of innate and B-cell-associated transcripts provide a useful example. Active tuberculosis is characterized by prominent interferon-inducible and myeloid blood programs [16,17], whereas the later appearance of B-cell-related expression during treatment was already observed by Cliff and colleagues [9]. In the present analysis, a reduction in estimated neutrophil representation accompanied a smaller rise in relative B-cell fraction. Because cell fractions are compositional, the apparent increase of one population can partly reflect contraction of another; neither deconvolution method measures absolute B-cell numbers. This does not make the B-cell signal biologically trivial. Independent cellular studies have described abnormalities in circulating B-cell phenotype and function during active tuberculosis and their improvement after treatment [51]. Rather, our findings show why evidence of B-cell-associated RNA recovery and evidence of B-cell numerical or functional reconstitution should be kept distinct. Both changes may occur, but bulk expression alone cannot establish their relative contributions.',
        'This distinction also refines the meaning of treatment-associated normalization. An early fall in inflammatory transcripts need not coincide with the slower restoration of adaptive immune features, and neither necessarily marks complete resolution of tissue disease. The staged transcriptional responses described by Bloom and Cliff [9,53] are consistent with different components of the host response evolving on different timescales. Thompson and colleagues further related blood RNA patterns to treatment outcome and pulmonary imaging, including evidence that inflammatory abnormalities may persist despite microbiological cure [10]. Taken together, these studies argue against treating a single blood score as a universal measure of bacterial clearance, immune competence and recovery of lung pathology. The value of a longitudinal host signature depends on which of these processes it is intended to measure.',
        'The consequences are especially important for inferred transcription-factor activity. Regulon-based methods estimate regulatory activity from the expression of target genes, rather than directly measuring transcription-factor action within a defined cell type [27,28]. When the relative abundance of target-expressing populations changes, a bulk regulon score may shift even without a corresponding change in intracellular regulation. In our data, treatment-associated PAX5, TBX21 and MAF effects were attenuated after accounting for estimated immune-cell composition. This interpretation is consistent with the distinct cellular localization of PAX5 in B-lineage cells and TBX21 in T/NK compartments observed in our single-cell analysis (Fig. 4), rather than providing evidence for a direct TBX21-PAX5 regulatory cascade. Importantly, attenuation after adjustment is a test of sensitivity to composition, not a causal estimate of how much recovery is mediated by cell abundance. Composition estimates and regulon scores are derived from overlapping expression information, and adjustment can remove genuine biological variation as well as confounding. The limited cross-cohort consistency of gene-specific attenuation further argues against assigning a stable cell-intrinsic fraction to individual genes.',
        'The same interpretive distinction applies to blood RNA biomarkers. A score can discriminate treatment states because it summarizes a broad change in immune-cell mixture, yet contribute little information about future outcome once that mixture is taken into account. The four-gene TBRI illustrates this problem: its close association with estimated B-cell fraction, including after accounting for treatment time, makes it a useful representation of the recovering blood expression state, but not an independently established measure of treatment success. The limited number of poor outcomes prevents a definitive assessment of incremental prognostic value, and the absence of a consistent early gain must not be interpreted as statistical equivalence. Importantly, this finding does not invalidate host RNA biomarkers in general. Thompson and colleagues demonstrated that blood RNA signatures provide information complementary to sputum-based Gene Xpert measurements and are associated with pulmonary inflammatory burden assessed by PET-CT [10]. Separately, transcriptional signatures have shown potential for predicting treatment outcomes in patients with tuberculosis and diabetes [18]. Signatures developed for diagnosis or future disease progression address different clinical questions [33,34]. A clinically useful treatment-monitoring marker therefore needs to demonstrate added value over readily available clinical, microbiological and cell-composition information, not merely track the passage of treatment time.',
        'There is a broader implication for studies of inflammation based on mixed-cell specimens. Longitudinal expression changes can arise from altered leukocyte trafficking, shifts in population abundance, changes in the activation state of existing cells, or combinations of these processes. These are different levels of explanation, and they imply different follow-up experiments. Composition-aware analysis can help reveal when an apparent regulatory signature is sensitive to the cellular mixture, but should not be used as an automatic procedure for removing every immune-cell-associated signal. If redistribution of leukocyte populations is itself part of the treatment response, indiscriminate adjustment risks subtracting the biology under investigation. The appropriate conclusion is therefore not that composition-independent signals are the only meaningful ones, but that the level of biological inference should match the type of measurement.',
        'Several features of the available data constrain that inference. The discovery contrasts share a baseline group, although directional testing in a separate cohort reduces reliance on this internally correlated comparison. The two deconvolution approaches operate on the same bulk transcriptomes and provide algorithmic sensitivity analyses rather than orthogonal cell counts; their estimates should also be interpreted cautiously across different assay platforms. The GSE89403 outcome group was small, and the upstream handling of deposited technical replicates could not be fully reconstructed from the processed matrix. All three longitudinal blood cohorts excluded participants with HIV, and clinical covariates were not consistently available. Finally, the lung and airway single-cell datasets provide cross-sectional cellular localization, not treatment-linked evidence of intracellular recovery. These limitations do not negate the reproducible treatment-associated expression pattern, but they restrict claims about its mechanism and clinical transportability.',
        'The next step is thus a more discriminating validation strategy rather than another search for a universal recovery gene or master regulator. Paired absolute leukocyte counts and flow-cytometric B-cell measurements, ideally combined with expression profiling of purified populations or donor-linked longitudinal single-cell data, would separate population turnover from within-cell transcriptional change. Larger prospective cohorts with serial microbiological measurements and, when feasible, imaging would then establish whether composition-sensitive and composition-independent RNA features convey complementary clinical information. By revisiting well-characterized treatment cohorts at these different levels, this study provides an interpretive framework for tuberculosis and related inflammatory conditions: transcriptional recovery is a robust observation, but its mechanistic and clinical meaning depends on distinguishing changes in the cells present from changes occurring within those cells.',
        'In conclusion, tuberculosis treatment-associated blood transcriptional recovery is reproducible, but its interpretation requires accounting for immune-cell composition, particularly when inferring transcription-factor activity and evaluating RNA-based biomarkers. The present evidence supports composition-associated recovery rather than absolute B-cell reconstitution, an intracellular causal mechanism, or established early clinical predictive utility.'
    ]
    # Locate the headings after the Results reordering and replace every paragraph
    # between Discussion and Funding, independent of fragile paragraph indices.
    current = list(document.paragraphs)
    discussion_heading = next(q for q in current if q.text.strip() == '4. Discussion')
    # The manuscript template stores front-matter metadata as a paragraph beginning
    # with "Author contributions" rather than a standalone Funding heading.
    funding_heading = next((q for q in current if q.text.strip() == 'Funding'), None)
    if funding_heading is None:
        funding_heading = next(q for q in current if q.text.strip().startswith('Author contributions'))
    node = discussion_heading._p.getnext()
    while node is not None and node is not funding_heading._p:
        nxt = node.getnext()
        node.getparent().remove(node)
        node = nxt
    for text in final_discussion:
        q = document.add_paragraph(text)
        q._p.getparent().remove(q._p)
        funding_heading._p.getparent().insert(funding_heading._p.getparent().index(funding_heading._p), q._p)
    # Keep the reproducibility statement concise and free of implementation details.
    for q in document.paragraphs:
        if 'processed result tables, and figure-generation code and derived result tables' in q.text:
            q.text = q.text.replace('processed result tables, and figure-generation code and derived result tables', 'derived result tables, and figure-generation code')
    for q in document.paragraphs:
        updated = remap(q.text)
        if updated != q.text: q.text = updated
    # The editorial pass may remove the original Funding block together with
    # the old Discussion tail. Restore the required funding statement directly
    # before Author contributions when it is absent.
    if not any(q.text.strip() == 'Funding' for q in document.paragraphs):
        author_heading = next((q for q in document.paragraphs
                               if q.text.strip() == 'Author contributions'), None)
        if author_heading is not None:
            h = document.add_paragraph('Funding', style='Heading 1')
            b = document.add_paragraph(
                'This work was supported by the China Postdoctoral Science Foundation '
                '(2020M671561) and the Jiangsu Postdoctoral Science Foundation (2020Z130).'
            )
            parent = author_heading._p.getparent()
            parent.insert(parent.index(author_heading._p), h._p)
            parent.insert(parent.index(author_heading._p), b._p)
    # Figure calls are remapped globally for the reordered Results sections.
    # This Discussion sentence refers specifically to the final localization
    # panel, which is Figure 4 in the submission legend block.
    for q in document.paragraphs:
        if 'distinct cellular localization of PAX5' in q.text:
            q.text = q.text.replace('(Fig. 5)', '(Fig. 4)')

def revise_supplement(document, figure_embed, add_body, add_heading):
    contents = next(p for p in document.paragraphs if p.text.startswith('Figure S4. Paired'))
    entry = contents.insert_paragraph_before('Figure S5. Cross-sectional airway single-cell context.', style='List Bullet')
    contents._p.addnext(entry._p)
    add_heading(document, 'Supplementary contextual results', 1)
    add_body(document, MOVED['network'])
    for text in MOVED['context']: add_body(document, text)
    add_heading(document, 'Figure S5', 1)
    for text in MOVED['bal']: add_body(document, text)
    add_body(document, MOVED['caption'])

def revise_cover(document):
    for p in document.paragraphs:
        if p.text.startswith('The central advance'):
            p.text = 'The study examines how immune-cell composition shapes the interpretation of reproducible blood transcriptional recovery. Independent gene-level directional testing supports the recovery pattern, while concordant relative-cell-fraction changes and attenuation of TF treatment effects after adjustment caution against interpreting bulk regulon recovery as intracellular regulation.'
        elif p.text.startswith('The analysis uses'):
            p.text = 'Three longitudinal cohorts have distinct discovery, score-development and outcome-association roles. The four-gene TBRI provides an example of the interpretive consequences: its association with estimated B-cell fraction persists after treatment-time adjustment, but evidence for incremental early outcome information is limited. Shared expression inputs and the absence of independent cell counts remain explicit limitations.'
        elif p.text.startswith('The manuscript is relevant'):
            p.text = 'The manuscript is relevant to inflammation research because it evaluates how changes in estimated leukocyte composition affect interpretation of inflammatory recovery, TF activity and RNA biomarkers. Cross-sectional single-cell observations and network priors provide context; detailed airway STAT4 results are confined to the supplement.'

