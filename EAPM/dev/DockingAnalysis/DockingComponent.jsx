import React, { useEffect, useState } from "react";

export default function DockingAnalysisView() {
  const {dockingFolder, modelsFolder} = window.parent.extensionData;

  console.log(window.parent.extensionData)

  const {href} = window.location

  return (
    <div>
      <p>Docking folder: {dockingFolder}</p>
      <p>Models folder: {modelsFolder}</p>
      <button
        onClick={() => {
          console.log(window.parent.extensionData)
        }}
      >
        Analyse
      </button>
    </div>
  );
}
