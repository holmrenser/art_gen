extern crate plotters;
use plotters::prelude::*;
use plotters::coord::cartesian::Cartesian2d;
use plotters::coord::types::RangedCoordi32;
use plotters::style::RGBColor;
use csv::ReaderBuilder;
use std::collections::HashMap;
use std::error::Error;
use serde::Deserialize;

type BitmapCartesianDrawingArea = DrawingArea<BitMapBackend<'static>, Cartesian2d<RangedCoordi32, RangedCoordi32>>;

const BACKGROUND: RGBColor = RGBColor(4,90,141);
const WHITE: RGBColor = RGBColor(255,247,251);
const BLUE: RGBColor = RGBColor(54,144,192);
const RED: RGBColor = RGBColor(227,26,28);
const YELLOW: RGBColor = RGBColor(255,237,160);

const SCALEFACTOR: i32 = 100;

fn bit_field(
    drawing_area: &BitmapCartesianDrawingArea
) -> () {
    for x in drawing_area.get_x_range() {
        for y in drawing_area.get_y_range() {
            match (x ^ y) % 19 {
                0|10 => {drawing_area.draw_pixel((x,y), &WHITE);()},
                3 => {drawing_area.draw_pixel((x,y), &YELLOW);()},
                1|11 => {drawing_area.draw_pixel((x,y), &BLUE);()},
                4|5 => {drawing_area.draw_pixel((x,y), &RED);()},
                _ => (),
            };
        }
    }
}

fn back(width: i32, height: i32) {
    let root = BitMapBackend::new("back.png", (width as u32, height as u32));
    let drawing_area = root
        .into_drawing_area()
        .apply_coord_spec(
            Cartesian2d::<RangedCoordi32, RangedCoordi32>
                ::new(-190i32..309i32, -438i32..271i32,(0..width, 0..height))
        );
    drawing_area.fill(&BACKGROUND);
    bit_field(&drawing_area);
}

#[derive(Debug, Deserialize)]
struct ScfSize {
    scf: String,
    len: i32,
}

#[derive(Debug, Deserialize)]
struct PafEntry {
    q_seqid: String,
    q_len: i32,
    q_start: i32,
    q_end: i32,
    strand: char,
    t_seqid: String,
    t_len: i32,
    t_start: i32,
    t_end: i32,
    n_match: i32,
    aln_len: i32,
    map_q: i32,
}

type ScaffoldSizes = HashMap<String, i32>;

fn get_scaf_sizes(path: &str) -> Result<ScaffoldSizes, Box<dyn Error>> {
    let mut sizes: ScaffoldSizes = HashMap::new();
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(path)?;
    
    for result in rdr.deserialize() {
        let record: ScfSize = result?;
        sizes.insert(record.scf, record.len);
    }
    Ok(sizes)
}

fn plot_alignments(
    drawing_area: &BitmapCartesianDrawingArea,
    pan_individual_sizes: ScaffoldSizes,
    tor_individual_sizes: ScaffoldSizes,
) -> Result<(), Box<dyn Error>> {
    let mut pan_cumulative_sizes: ScaffoldSizes = HashMap::new();
    let mut tor_cumulative_sizes: ScaffoldSizes = HashMap::new();

    drawing_area.fill(&YELLOW);

    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_path("./alignments/pantor.minimap.sorted.paf")?;

    let mut pan_cumulative_offset = 0;
    let mut tor_cumulative_offset = 0;
    for result in rdr.deserialize() {
        let aln: PafEntry = result?;
        if aln.aln_len < 1000 {
            continue
        }
        let pan_offset = match pan_cumulative_sizes.get(&aln.t_seqid) {
            Some(offset) => *offset,
            None => {
                let pan_scf_size = *pan_individual_sizes.get(&aln.t_seqid)
                    .ok_or("KEY NOT FOUND")?;
                pan_cumulative_sizes.insert(aln.t_seqid,
                    pan_cumulative_offset + pan_scf_size);
                pan_cumulative_offset = pan_cumulative_offset + pan_scf_size;
                pan_cumulative_offset - pan_scf_size
            },
        };
        let tor_offset = match tor_cumulative_sizes.get(&aln.q_seqid) {
            Some(offset) => *offset,
            None => {
                let tor_scf_size = *tor_individual_sizes.get(&aln.q_seqid)
                    .ok_or("KEY NOT FOUND")?;
                tor_cumulative_sizes.insert(aln.q_seqid,
                    tor_cumulative_offset + tor_scf_size);
                tor_cumulative_offset = tor_cumulative_offset + tor_scf_size;
                tor_cumulative_offset - tor_scf_size
            },
        };
        
        let colour = match aln.strand {
            '+' => &RED,
            '-' => &BLUE,
            _ => &YELLOW
        };
        for i in 0..(aln.aln_len) {
            let (x,y) = (
                (pan_offset + aln.t_start + i) / SCALEFACTOR,
                (tor_offset + aln.q_start + i) / SCALEFACTOR
            );
            drawing_area.draw_pixel((x, y), colour);
        };
    };
    Ok(())
}

fn front(width: i32, height: i32) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new("front.png", (width as u32, height as u32));
    let margin_top = 40;
    let margin_bottom = 40;
    let margin_left = 40;
    let margin_right = 40;
    let center_width = width - margin_left - margin_right;
    let center_height = height - margin_top - margin_bottom;
    let main_drawing_area = root.into_drawing_area();
    main_drawing_area.fill(&BACKGROUND);
    let areas = main_drawing_area.split_by_breakpoints(
        [margin_left, width - margin_right],
        [margin_top, height - margin_bottom]
    );

    bit_field(&main_drawing_area.apply_coord_spec(
        Cartesian2d::<RangedCoordi32, RangedCoordi32>
            ::new(-width..0, 0..height,(0..width, 0..height))
    ));

    let (mid1, mid2) = areas[4].split_vertically(
        center_height - center_width
    );
    let title_areas = mid1.split_by_breakpoints(
        [20, center_width - 20],
        [20, center_height - center_width - 40]
    );
    title_areas[4].fill(&BACKGROUND);

    let pan_individual_sizes = get_scaf_sizes(
        &"./alignments/PanWU01x14_asm01.scf.masked.duprm.scaffold_lengths.tsv"
    )?;
    
    let tor_individual_sizes = get_scaf_sizes(
        &"./alignments/TorRG33x02_asm01.scf.masked.duprm.scaffold_lengths.tsv"
    )?;
    
    let alignment_area = mid2.apply_coord_spec(
        Cartesian2d::<RangedCoordi32, RangedCoordi32>
            ::new(
                0..(pan_individual_sizes.values().sum::<i32>() / SCALEFACTOR),
                0..(tor_individual_sizes.values().sum::<i32>() / SCALEFACTOR),
                (
                    margin_left + 10..(width - margin_right - 10),
                    (height - margin_bottom - center_width + 10)..(height - margin_bottom - 10)
                )
            )
    );
    alignment_area.fill(&YELLOW);
    plot_alignments(&alignment_area, pan_individual_sizes, tor_individual_sizes)
}

fn main() {
    let width: i32 = 499;
    let height: i32 = 709;
    back(width, height);
    front(width, height);
}
