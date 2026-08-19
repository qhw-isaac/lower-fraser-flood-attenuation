# Walkthrough clips

`index.html` shows these six clips to phone visitors, in this order. Each is
fetched only as it nears the viewport, and tapping one opens it full screen.
A file that is absent is replaced at runtime by a labelled placeholder carrying
the clip's description, so the page stays readable while a recording is
outstanding.

| File | Section | What it shows |
|------|---------|---------------|
| `community-select.mp4` | 1 Community view | Abbotsford is picked, the view pulls back, the channels draining toward it are drawn in |
| `community-compare.mp4` | 1 Community view | Chilliwack is added as the second community, the ground upstream of both turns gold |
| `upstream-ecosystems.mp4` | 2 Upstream ecosystems | A watershed is selected, then the shading is switched across the four readings |
| `supply-select.mp4` | 3 Supply & demand | A watershed is selected on the supply side, the chain carrying its water to the outlet lights up |
| `demand-select.mp4` | 3 Supply & demand | Neighbourhoods are selected on the demand side, each reporting its exposed area |
| `weights.mp4` | 3 Supply & demand | The exposure weight sliders are dragged and the ranking re-shades |

Each clip carries a poster of the same basename with a `.webp` extension, loaded
at the same time as the clip.

## Preparing a recording

Masters are kept in `context/recordings/`, which is outside the published site
and outside git. Screen recordings arrive as 2880 x 1536 QuickTime `.mov` at
60 fps. Browsers make no commitment to the QuickTime container, so convert
before publishing:

```bash
ffmpeg -i master.mov -an -vf "fps=30,scale=1280:-2" \
  -c:v libx264 -profile:v main -pix_fmt yuv420p -crf 26 \
  -movflags +faststart clip.mp4

ffmpeg -ss 0.2 -i master.mov -frames:v 1 -vf "scale=1280:-2" clip.webp
```

- `libx264` + `yuv420p` is the one combination every browser decodes. HEVC,
  which newer Macs record by default, plays in Safari and not in Chrome.
- `+faststart` moves the index to the front so playback begins before the file
  finishes downloading.
- `-an` drops audio, which the page never plays.

Any aspect ratio works. The page reads each clip's own dimensions once its
metadata loads, and reserves 15:8 until then, so a set that is not 15:8 should
have that default changed to match.

Preview the page from a desktop at `index.html?guide`.
