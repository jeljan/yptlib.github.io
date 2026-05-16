import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

const hover = JSON.parse(
  await readFile(new URL("../assets/data/os/hover/Z1688207185.json", import.meta.url), "utf8")
).rows;

assert.deepEqual(
  hover["4316"],
  [
    [87.572741, 90.868031, 80.57505, 67.172039],
    [98.842806, 113.579434],
  ],
  "Q13263 Y369 hover should use the DMSO channels and Z1688207185 channels from plex 27"
);

assert.deepEqual(
  hover["13333"],
  [
    [77.119533, 73.355951, 32.101671, 50.970232],
    [64.918589, 73.517053],
  ],
  "Q5K4L6 Y673 hover should use the DMSO channels and Z1688207185 channels from plex 27"
);
