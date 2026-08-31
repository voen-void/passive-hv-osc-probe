# Passive High Voltage oscilloscope probe

> [!IMPORTANT]
> **UNTESTED!**
> **This project is still under construction and has not yet been electrically tested. The stated voltage ratings are design targets, not verified operating limits. It is not recommended to build this currently!**

![A drawing of the schematic of the probe](/images/schematic-basic.jpeg)

[Schematic pdf](/images/hv-probe-passive-schematic.pdf) included in the repo.

*Summary*
A 1000:1 passive high-voltage oscilloscope probe designed for measurements up to 20 kV peak. The probe uses a compensated resistive/capacitive divider, a shielded and resin-potted PCB construction, and is designed to fit inside a 37.6 mm ID metal tube (1-1/2" tubing will do fine here).

## TLDR
| Parameter                    | Value             |
| ---------------------------- | ----------------- |
| Attenuation                  | 1000:1            |
| Maximum recommended voltage  | 20 kV peak        |
| Divider                      | 50 MΩ / 50 kΩ     |
| Probe body                   | Resin-potted PCB  |
| Shield                       | Brass/copper tube |
| Shield inner diameter        | 37.6 mm           |
| PCB                          | 2-layer           |
| Recommended potting resin εr | ≤ 6               |
| Compensation                 | Adjustable        |
| Status                       | **Not yet tested**|

The design is based on a [EEVBlog Video](https://www.youtube.com/watch?v=jUvSP3BQpvs). The electrical design is done in [KiCad](https://www.kicad.org/) 10.0.5, and the mechanical design is done in [onshape](https://www.onshape.com/en/). The total cost of a single probe adds up to around **100 USD**, though it will be larger because not all parts can be bought in small enough quantities (eg. the resin, the PCBs \[usually 5 minimum\] and the shielding tube).

![A rendering of the front of the PCB of the probe](/images/pcb-render-front.jpeg)
![A rendering of the back of the PCB of the probe](/images/pcb-render-back.jpeg)

> [!NOTE]
> The PCBs for this project have been provided by [PCBWay](https://www.pcbway.com/).

## What is it?
This project tackles the design of a simple high voltage oscilloscope probe based on the principles of a voltage divider. The simplest schematic would look like this:


<!--
![A drawing of a simple resistive divider](/images/divider.jpg)
-->

<img src="/images/divider.jpg" width="30%" />


It's a bit more complex than that if we want to probe AC signals of considerable frequency, as the parasitic capacitances need to be compensated to achieve better bandwidth. For this reason the resistors forming the voltage divider are bypassed by capacitors which increase the AC transmittance of the circuit. It should be noted that for desired probe behavior the transmittance of just the resistors should match that of just the capacitors. If we cared only to measure DC, then the capacitors could be ommited entirely.

At the low voltage end two footprints for capacitors are provided which allow simple compensation as you will need to adjust that capacitance to achieve a flat frequency response. Capacitor values are also proposed, but the final capacitor value may differ from the proposed value because the parasitic capacitance depends on the PCB, components, potting compound and mechanical construction.
For adjusting the DC and low frequency transmittance a 10 turn potentiometer is provided. It should be set to 1 M (full range) by default. While calibrating adjust it to reach desired 1000:1 division ratio. Then the capacitor should be adjusted to achieve flat frequency response. A proposed calibration procedure follows;

### Calibration

The proposed calibration procedure:
1. Set the 10-turn potentiometer to 1 MΩ.
2. Apply a suitable low-voltage signal.
3. Adjust the potentiometer until the DC attenuation is 1000:1.
4. Apply a square-wave calibration signal.
5. Adjust the compensation capacitor(s) for the flattest frequency response.
6. Verify the response at multiple frequencies.

## Safety

> [!CAUTION]
> High voltages can easily be lethal. Only deal with them if you are qualified. You are building and using this project at your own risk! This is experimental high-voltage equipment.

1. Do not use the probe on energized equipment unless you are qualified to work with high voltage.
2. The 20 kV value is currently a design target, not a verified rating.
3. Never rely on the probe alone to provide isolation.
4. Use an appropriate current-limited HV source during testing.
5. Ensure the oscilloscope and all connected equipment are rated appropriately.
6. The metal shield must be properly grounded.
7. Inspect the probe for mechanical damage before use.
8. Do not use a probe that has experienced an internal arc-over.
7. Potting compound can make internal failures difficult or impossible to inspect.

The maximum recommended voltage which this probe can measure is 20 kV **peak**. The components should handle 50 kV, but the final breakdown voltage depends on the quality of your craftsmanship - primarly the resin used, how clean the elements were when casting them in resin, if there are any air bubbles trapped inside the resin (**there shouldn't be any!**). A safety margin must be used!

### Internal arc-over protection

For increased safety a microwave fuse may be added in series with the input of the probe. This way if an over-current condition occurs - indicating an internal arc over - the fuse will open indicating thus. Microwave fuses are the highest voltage rated fuses I could cheaply find, and should do the job. A lower current fuse would be ideal, but those only come with voltage ratings up to 250 VAC. The reason for that fuse is to protect the probe in the case of an internal arc-over. If an arc over was to happen inside the probe, excessive current would be drawn from DUT (if DUT can provide it). Sustained internal arcing might lead to pressure buildup inside the probe body, which due to it's construction could become very dangerous. In an extreme scenario the probe could explode causing harm to everything around it. For this reason it is recommended to use some sort of fuse in series with the probe input, and it it blows, the probe should no longer be used and should be marked as faulty and discarded. High voltage apparatus potted in hard resin is known to pose this kind of danger and should be handled accordingly.


## Physical construction

![A drawing with the dimensions of the PCB of the probe](/images/pcb-basic.jpeg)

[A drawing with the dimensions](/images/hv-probe-passive-dimensions.pdf) of the PCB included in the repo.

The probe requires some special attention when it comes to construction. Since the divider consists of such high value resistors (50 MOhm to 50 kOhm divider), it's sensitive to interference. For this reason the probe is designed to fit in a shield in the form of a brass or copper tube of the inner diameter of 37,6 mm (1.48 in). This shield shall be grounded to protect the sensitive signal inside. To increase dielectric breakdown a special resin must be used, and the PCB inside the shield should be cast in it. Since the voltages measured are significant, leaving the PCB in the free air wont cut it. Try to pick a resing with as low dielectric constant as possible. [HUNTSMAN CW5620 BLUE](https://products.huntsman.com/products/arathane-cw-5620-hy-5610) was used, which has a relative dielectric constant of 6. It is recommended to not use resins of higher relative dielectric constant as it increase the parasitic capacitances, which will degrade the bandwidth of the probe cutting down high frequency transmittance.

### Cost of building a single probe

| Element (MPN)                | Reference     | Qty.  | Price \[USD\] |
| --------------------------   | ------------- | ----- | ----------- |
| PCB                          | -             | 1     | 4.2         |
| 564RC0GAJ602EJ240J           | C1-C5, C8-C12 | 10    | 30.27       |
| VJ1206A222FFBAT              | C6, C7        | 2     | 2.55        |
| 112413                       | J2            | 1     | 8.22        |
| VR68000001005JAC00           | R1-R5         | 5     | 7.02        |
| CR1206-FX-5602ELF            | R6            | 1     | 0.1         |
| 3296Z-1-105LF                | RV1           | 1     | 2.45        |
| HUNTSMAN CW5620 BLUE Resin   | -             | 732 g | 25          |
| Copper shielding tube        | -             | 15 cm | 17          |
| 3DP filament cost            | -             | -     | 2           |

Total: 98.81 USD without shipping for a single probe. Both the PCBs, the resin and the shielding tube only come in larger quantities than needed for the project, thus the actual price may be higher depending on what you are able to source.

### PCB

The PCB has been designed to be a simple two-layer board. Soldermask has been removed from the part of the PCB that shall be submerged in isolating resin. The part of the PCB that has soldermask left should be left sticking out of the resin to allow for adjustments of the compensation capacitors and resistor if needed. A solder pad connected to probe ground is exposed to facilitate an easy connection to the shield.

The PCB has been designed with a specific kind of HV terminal construction. A large portion of copper is exposed at the HV end of the probe with a long distance separating it from the shield, to increase creepage distance. This part of the probe (tapered end) should be exposed from the shield and should be cast in resin, preferably in a ribbed shape for further increase of creepage distance. Before casting the PCB a brass piece shall be machined such that it fits over the straight tip of the tapored end of the PCB and presents a female thread at the end, into which a screw may be threaded into to make electrical contact. This piece should be soldered to the tapored tip of the PCB before casting. It should be noted that **no excessive force** should be used on this part of the probe in relation to the probe body. To make fastening easier the exposed part of the brass piece may be shaped in a way which enables holding it with an appropirately sized wrench. In such case a wrench may be used to hold the HV terminal of the probe down while fastening an attachment to it. The end piece should only be tightened finger-tight.

#### The manufactured PCB

I had the PCB manufactured by [PCBWay](https://www.pcbway.com/). The process was very simple: I uploaded the gerbers, placed a Quick Order and the PCBs (five of them!) arrived at my doorstep in a couple of days. I was pleasantly surprised by the packing quality and the PCB manufacture quality was excellent, picture below!

![A photo of the manufactured PCB](/images/pcb-made.jpg)

After receiving the PCB i shortly proceeded to populating it with components. [Assembly drawings](/images/hv-probe-passive-assembly.pdf) included in the repo.

### Notes on elements used

![An assembly drawing of the PCB of the probe](/images/pcb-assembly-front.jpeg)

The resistors used for the high voltage part of the voltage divider are special high voltage resistors, which may operate up to 10 kV across them. The capacitors used in the divider are all C0G / NP0 dielectric capacitors, for increased  temperature stability. other dielectrics may be used, but it's strongly recommended to use Class I dielectric materials, preferably C0G / NP0 as it will cause the smallest drift in alignment over a wide range of temperatures.

### Casting in resin

To cast the probe in resin a CAD model for a 3D-printable mold has been made in [onshape](https://cad.onshape.com/documents/634b2c7804eb6959816cbee4/w/1f30c7d8f223e5e8fad3fc6d/e/e251d38d61b503254921d7bb). 

<!--
![A rendering of one half of the 3D-printable mold](/images/mold-side-view.png)
-->

<img src="/images/mold-side-view.png" width="50%" />

Recommended material for the mold is PETG, as it's recommended for molding and should allow for easier separation after casting. The two mold halves are to be screwed together by 18 70 mm long M4 bolts. Wax was used as mold release. The resin should be poured down into the metal tube, with the tip of the probe facing downwards. Resin should be added until it fills up the mold up to the top 3 mm venting holes.