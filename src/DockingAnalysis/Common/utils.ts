export async function postData(url = "", data = {}) {
  const response = await fetch(url, {
    method: "POST", // *GET, POST, PUT, DELETE, etc.
    headers: {
      "Content-Type": "application/json",
      Accept: "application/json",
    },
    body: JSON.stringify(data), // body data type must match "Content-Type" header
  });
  return response.json(); // parses JSON response into native JavaScript objects
}

export async function postOctet(url = "", data = {}) {
  const response = await fetch(url, {
    method: "POST", // *GET, POST, PUT, DELETE, etc.
    headers: {
      "Content-Type": "application/octet-stream",
      Accept: "application/octet-stream",
    },
    body: JSON.stringify(data), // body data type must match "Content-Type" header
  });
  return response.json(); // parses JSON response into native JavaScript objects
}

// Example POST method implementation:
export async function getData(url = "") {
  // Default options are marked with *
  const response = await fetch(url, {
    method: "GET", // *GET, POST, PUT, DELETE, etc.
  });
  return response.json(); // parses JSON response into native JavaScript objects
}

export function getFileForMolstar(data: {
  contents: string;
  filename: string;
}): File {
  // Use the parent Blob and File definitions. If not,
  // molstar will fail to load the file because wont detect
  // that the Blob is instance of Blob
  // window.Blob != parent.Blob
  // @ts-ignore
  const blob = new parent.Blob([data.contents]);
  // @ts-ignore
  const file = new parent.File([blob], data.filename ?? "Molecule.pdb");

  return file;
}
