pub mod haplotree;
pub mod haplotype;

use itertools::Itertools;

use crate::utils::seq::SeqRegion;
use crate::variant::Var;

// fn phase_variants(variants:&Vec<Var>) {
//     for (tid, group) in &variants.into_iter().group_by(|&var| var.tid) {
//         let mut positions = group.map(|var| var.pos).collect_vec();
//         positions.sort_unstable();
//         phase_region(tid, &positions);
//     }
// }


fn phase_region(seq_region:&SeqRegion) {
    todo!()
}

    // def separate_references(self):
    //     for reference_id in self.reference_lengths:
    //         reference_positions = self.variant_positions[reference_id]
    //         reference_phasesets = dict()
    //         reference_variant_positions = set()
    //         for reference_start, reference_end in self.reference_intervals[reference_id]:
    //             first = bisect.bisect_left(reference_positions, reference_start)
    //             last = bisect.bisect_right(reference_positions, reference_end-1)
    //             region_phasesets, region_positions = self._separate_reference_region(reference_id, reference_start, reference_end, reference_positions[first:last])
    //             reference_phasesets |= region_phasesets
    //             reference_variant_positions |= region_positions
    //         self.phasesets[reference_id] = reference_phasesets
    //         self.variant_positions[reference_id] = reference_variant_positions