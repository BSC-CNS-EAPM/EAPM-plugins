import React from "react";
import ReactDOM from "react-dom/client";
import DockingAnalysisView from "./DockingComponent";

// Include the default CSS
import "./docking_analysis.css";

// Create the container
const container = document.createElement("div");
container.id = "root-page";
document.body.appendChild(container);

ReactDOM.createRoot(container).render(
  <div className="root-page">
    <DockingAnalysisView />
  </div>
);
