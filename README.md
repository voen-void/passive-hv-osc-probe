# Passive High Voltage oscilloscope probe
KiCad design files for a classic passive 1000:1 high voltage oscilloscope probe.

Design based on an old [EEVBlog Video](https://www.youtube.com/watch?v=jUvSP3BQpvs).

As in that video, the probe has been designed to fit in a shielding tube. The probe has been designed to be cast in resin for increased dielectric breakdown strength.

![A drawing of the schematic of the probe](/images/schematic-basic.jpeg)

![A rendering of the front of the PCB of the probe](/images/pcb-render-front.jpeg)
![A rendering of the back of the PCB of the probe](/images/pcb-render-back.jpeg)

> [!NOTE]
> **This project is still under construction! It is not recommended to build it yet, as it has not yet been tested.**
> *The PCBs for this project have been provided by PCBWay.*

## What is it?
This project tackles the design of a simple high voltage oscilloscope probe based on the principles of a voltage divider. Its a bit more complex than that, as the parasitic capacitances need to be compensated to achieve better bandwidth. For this reason the resistors forming the voltage divider are bypassed by capacitors which increase the AC transmittance of the circuit. It should be noted that for desired probe behavior the tranmittance of just the resistors should match that of just the capacitors. 

At the low voltage end two footprints for capacitors are provided which allow simple compensation as you will need to adjust that capacitance to achieve a flat frequency response. Capacitor values are also proposed, but might need to be changed.
For adjusting the DC and low frequency transmittance a 10 turn potentiometer is provided. It should be set to 1 M (full range) by default. While calibrating adjust it to reach desired 1000:1 division ratio. Then the capacitor should be adjusted to achieve flat frequency response.


## TLDR What you need to know if you want to build it
1. Reading the whole document is strongly recommended.
2. The maximum recommended voltage which this probe can measure is 20 kV PEAK. The components should handle 50 kV, but the final breakdown voltage depends on the quality of your craftsmanship - primarly the resin used, how clean the elements were when casting them in resin, if there are any air bubbles trapped inside the resin (there shouldnt be any). A safety marigin must be used!
3. The PCB has been designed to fit inside of a brass tube with the inner diameter of 37,6 mm (1.48 in)
4. The total cost of the project is around XX $.
5. Use a resin with a relative dielectric constant no higher than 6.

> [!CAUTION]
> High voltages can easily be lethal. Only deal with them if you are qualified. You are building and using this project at your own risk!

## Physical construction

![A drawing with the dimensions of the PCB of the probe](/images/pcb-basic.jpeg)

The probe requires some special attention when it comes to construction. Since the divider consists of such high value resistors (50 MOhm to 50 kOhm divier), its sensitive to interference. For this reason the probe is designed to fit in a shield in the form of a brass or copper tube of the inner diameter of 37,6 mm (1.48 in). This shield shall be grounded to protect the sensitive signal inside. To increase dielectric breakdown a special resin must be used, and the PCB inside the shield should be cast in it. Since the voltages measured are significant, leaving the PCB in the free air wont cut it. Try to pick a resing with as low dielectric constant as possible. HUNTSMAN CW5620 BLUE was used, which has a relative dielectric constant of 6. It is recommended to not use resins of higher relative dielectric constant as it increase the parasitic capacitances, which will degrade the bandwidth of the probe cutting down high frequency transmittance.

### PCB

The PCB has been designed to be a simple two-layer board. Soldermask has been removed from the part of the pcb that shall be submerged in isolating resin. The part of the PCB that has soldermask left should be left sticking out of the resin to allow for adjustments of the compensation capacitors and resistor if needed. A solder pad connected to probe ground is exposed to facilitate an easy connection to the shield.

The PCB has been designed with a specific kind of HV terminal construction. A large portion of copper is exposed at the HV end of the probe with a long distance separating it from the shield, to increase creepage distance. This part of the probe (tapored end) should be exposed from the shield and should be cast in resin, preferably in a ribbed shape for further increase of creepage distance. Before casting the PCB a brass piece shall be machined such that it fits over the straight tip of the tapored end of the PCB and presents a female thread at the end, into which a screw may be threaded into to make electrical contact. This piece should be soldered to the tapored tip of the PCB before casting. It should be noted that NO EXCESSIVE FORCE should be used on this part of the probe in relation to the probe body. To make fastening easier the exposed part of the brass piece may be shaped in a way which enables holding it with an appropirately sized wrench. In such case a wrench may be used to hold the HV terminal of the probe down while fastening an attachment to it. The end piece should only be tightened finger-tight.

### Notes on elements used

The resistors used for the high voltage part of the voltage divider are special high voltage resistors, wchich may operate up to 10 kV across them. The capacitors used in the divider are all C0G / NP0 dielectric capacitors, for increased  temperature stability. other dielectrics may be used, but its strongly recommended to use Class I dielectric materials, preferably C0G / NP0 as it will cause the smallest drift in alignment over a wide range of temperatures.

### Notes on internal arc-overs

For increased safety a microwave fuse may be added in series with the input of the probe. Microwave fuses are the highest voltage rated fuses I could cheaply find, and should do the job. A lower current fuse would be ideal, but those only come with voltage ratings up to 250 VAC. The reason for that fuse is to protect the probe in the case of an internal arc - over. If an arc over was to happen inside the probe, excessive current would be drawn from DUT (if DUT can provide it). Sustained internal arcing might lead to pressure buildup inside the probe body, which due to its construction could become very dangerous. In an extreme scenario the probe could violently explode causing much harm to everything around it. For this reason it is recommended to use some sort of fuse in series with the probe input, and it it blows, the probe should no longer be used and should be marked as faulty and discarded. High voltage apparatus potted in hard resin is known to pose this kind of danger and should be handled accordingly.
